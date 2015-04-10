
#include <math.h>
#include <stdio.h>
#include <time.h>
#include "cow.h"


/*
 * Macro for a three-dimensional loop over all interior cells
 * =====================================================================
 */
#define FOR_ALL_INTERIOR(N1, N2, N3)		\
  for (int i=2; i<N1+2; ++i)			\
    for (int j=2; j<N2+2; ++j)			\
      for (int k=2; k<N3+2; ++k)


/*
 * Macro to calculate linear index of (i,j,k,m) ... m goes from 1, not 0
 * =====================================================================
 */
#define IND(i,j,k,m) ((i)*si + (j)*sj + (k)*sk) + (m-1)


/*
 * Macros to calculate finite differences
 * =====================================================================
 */

#define GRADC2(F,s) ((+1*(F)[+1*s] +		\
		      -1*(F)[-1*s]) / 2.0)

#define GRADC4(F,s) ((+1*(F)[-2*s] +		\
		      -8*(F)[-1*s] +		\
		      +8*(F)[+1*s] +		\
		      -1*(F)[+2*s]) / 12.0)

#define CROSS(E,B) {0.0,				\
      (E)[2]*(B)[3]-(E)[3]*(B)[2],			\
      (E)[3]*(B)[1]-(E)[1]*(B)[3],			\
      (E)[1]*(B)[2]-(E)[2]*(B)[1]}			\

#define DOT(E,B) ((E)[1]*(B)[1] + (E)[2]*(B)[2] + (E)[3]*(B)[3])

#define MAX3(a,b,c) (a>b && b>c ? a : (b > c ? b : c))
#define MIN3(a,b,c) (a<b && b<c ? a : (b < c ? b : c))


enum FfeSimParameter {
  FFE_OHMS_LAW_VACUUM,
  FFE_OHMS_LAW_FORCE_FREE
} ;


struct ffe_measure
{
  double total_electric_energy;
  double total_magnetic_energy;
  double total_magnetic_helicity;
} ;


struct ffe_status
{
  int iteration;
  int checkpoint_number;
  double time_simulation;
  double time_final;
  double time_step;
  double kzps;
} ;


struct ffe_sim
{
  cow_domain *domain;
  cow_dfield *electric[8]; /* 0-4: E, 5-8: dtE */
  cow_dfield *magnetic[8];

  int Ni, Nj, Nk;
  double grid_spacing[4];
  double cfl_parameter;

  enum FfeSimParameter flag_ohms_law;

  struct ffe_status status;
  struct ffe_measure measure;
} ;



/*
 * Initialze a new simulation instance
 * =====================================================================
 */
void ffe_sim_init(struct ffe_sim *sim)
{
  sim->domain = cow_domain_new();

  int Ni = 8;
  int Nj = 8;
  int Nk = 64;

  sim->Ni = Ni;
  sim->Nj = Nj;
  sim->Nk = Nk;
  sim->cfl_parameter = 0.1;
  sim->flag_ohms_law = FFE_OHMS_LAW_VACUUM;
  //sim->flag_ohms_law = FFE_OHMS_LAW_FORCE_FREE;

  sim->status.iteration = 0;
  sim->status.checkpoint_number = 0;
  sim->status.time_simulation = 0.0;
  sim->status.time_step = 0.0;


  sim->grid_spacing[0] = 0.0;
  sim->grid_spacing[1] = 1.0 / Ni;
  sim->grid_spacing[2] = 1.0 / Nj;
  sim->grid_spacing[3] = 1.0 / Nk;

  cow_domain_setndim(sim->domain, 3);
  cow_domain_setsize(sim->domain, 0, Ni);
  cow_domain_setsize(sim->domain, 1, Nj);
  cow_domain_setsize(sim->domain, 2, Nk);
  cow_domain_setguard(sim->domain, 2);
  cow_domain_commit(sim->domain);

  for (int n=0; n<8; ++n) {
    sim->electric[n] = cow_dfield_new();
    sim->magnetic[n] = cow_dfield_new();

    cow_dfield_setname(sim->electric[n], "electric");
    cow_dfield_setname(sim->magnetic[n], "magnetic");

    cow_dfield_setdomain(sim->electric[n], sim->domain);
    cow_dfield_setdomain(sim->magnetic[n], sim->domain);

    cow_dfield_addmember(sim->electric[n], "E1");
    cow_dfield_addmember(sim->electric[n], "E2");
    cow_dfield_addmember(sim->electric[n], "E3");
    cow_dfield_addmember(sim->magnetic[n], "B1");
    cow_dfield_addmember(sim->magnetic[n], "B2");
    cow_dfield_addmember(sim->magnetic[n], "B3");

    cow_dfield_commit(sim->electric[n]);
    cow_dfield_commit(sim->magnetic[n]);
  }
}



/*
 * Finalize a simulation instance
 * =====================================================================
 */
void ffe_sim_free(struct ffe_sim *sim)
{

  for (int n=0; n<8; ++n) {
    cow_dfield_del(sim->electric[n]);
    cow_dfield_del(sim->magnetic[n]);
  }

  cow_domain_del(sim->domain);
}



/*
 * Evaluate initial data
 * =====================================================================
 */
void ffe_sim_initial_data(struct ffe_sim *sim)
{
  int Ni = cow_domain_getnumlocalzonesinterior(sim->domain, 0);
  int Nj = cow_domain_getnumlocalzonesinterior(sim->domain, 1);
  int Nk = cow_domain_getnumlocalzonesinterior(sim->domain, 2);
  int k0 = cow_domain_getglobalstartindex(sim->domain, 2);
  int si = cow_dfield_getstride(sim->electric[0], 0);
  int sj = cow_dfield_getstride(sim->electric[0], 1);
  int sk = cow_dfield_getstride(sim->electric[0], 2);
  double *E = cow_dfield_getdatabuffer(sim->electric[0]);
  double *B = cow_dfield_getdatabuffer(sim->magnetic[0]);
  double dz = sim->grid_spacing[3];

  FOR_ALL_INTERIOR(Ni, Nj, Nk) {

    double z = dz * (k + k0 - 2);

    int m = IND(i,j,k,0);

    E[m+1] = sin(2 * M_PI * z);
    E[m+2] = 0.0;
    E[m+3] = 0.0;

    B[m+1] = 0.0;
    B[m+2] = sin(2 * M_PI * z);
    B[m+3] = 1.0;

  }

  cow_dfield_syncguard(sim->electric[0]);
  cow_dfield_syncguard(sim->magnetic[0]);
}



/*
 * Evalute the current based on E, B, and gradients
 * =====================================================================
 */
void ffe_sim_ohms_law(struct ffe_sim *sim,
		      double E[4], double rotE[4], double divE,
		      double B[4], double rotB[4], double divB, double J[4])
{
  switch (sim->flag_ohms_law) {

  case FFE_OHMS_LAW_VACUUM:
    J[1] = 0.0;
    J[2] = 0.0;
    J[3] = 0.0;
    break;
  case FFE_OHMS_LAW_FORCE_FREE: /* eqn 11: Pfeiffer (2013) */
    {
      double S[4] = CROSS(E, B);
      double B2 = DOT(B, B);
      double Jb = DOT(B, rotB) - DOT(E, rotE);
      double Js = divE;

      J[1] = (B[1] * Jb + S[1] * Js) / B2;
      J[2] = (B[2] * Jb + S[2] * Js) / B2;
      J[3] = (B[3] * Jb + S[3] * Js) / B2;
    }
    break;
  }
}



/*
 * Advance the simulation by one Runge-Kutta substep
 * =====================================================================
 */
void ffe_sim_advance_rk(struct ffe_sim *sim, int RKstep)
{
#define D1(F,c) GRADC2(F+m+c,si)/dx
#define D2(F,c) GRADC2(F+m+c,sj)/dy
#define D3(F,c) GRADC2(F+m+c,sk)/dz

  double RKparam_array[4] = {0.0, 0.5, 0.5, 1.0};
  double RKparam = RKparam_array[RKstep];

  int Ni = cow_domain_getnumlocalzonesinterior(sim->domain, 0);
  int Nj = cow_domain_getnumlocalzonesinterior(sim->domain, 1);
  int Nk = cow_domain_getnumlocalzonesinterior(sim->domain, 2);

  int si = cow_dfield_getstride(sim->electric[0], 0);
  int sj = cow_dfield_getstride(sim->electric[0], 1);
  int sk = cow_dfield_getstride(sim->electric[0], 2);

  double *E0 = cow_dfield_getdatabuffer(sim->electric[0]);
  double *B0 = cow_dfield_getdatabuffer(sim->magnetic[0]);
  double *E  = cow_dfield_getdatabuffer(sim->electric[RKstep]); /* WRONG!!! */
  double *B  = cow_dfield_getdatabuffer(sim->magnetic[RKstep]);

  double *dtE = cow_dfield_getdatabuffer(sim->electric[RKstep+4]);
  double *dtB = cow_dfield_getdatabuffer(sim->magnetic[RKstep+4]);

  double dx = sim->grid_spacing[1];
  double dy = sim->grid_spacing[2];
  double dz = sim->grid_spacing[3];
  double dt = sim->status.time_step * RKparam;

  /* ===========================================================================
   * Fill in the n-th (n=0,1,2,3) time-derivative register, reading from the
   * n-th field register
   * ===========================================================================
   */
  FOR_ALL_INTERIOR(Ni, Nj, Nk) {

    int m = IND(i,j,k,0);

    double divE = D1(E,1) + D2(E,2) + D3(E,3);
    double divB = D1(B,1) + D2(B,2) + D3(B,3);    
    double rotE[4] = {0, D2(E,3) - D3(E,2), D3(E,1) - D1(E,3), D1(E,2) - D2(E,1)};
    double rotB[4] = {0, D2(B,3) - D3(B,2), D3(B,1) - D1(B,3), D1(B,2) - D2(B,1)};
    double J[4];


    /* Hyperbolicity terms, eqn 48-49: Pfeiffer (2013) */
    double B2 = DOT(&B[m], &B[m]);
    double gradEdotB[4] = {0, 0, 0, 0}; /* TODO */

    double S[4]   = CROSS(&E[m], &B[m]);
    double ct1[4] = CROSS(&E[m], gradEdotB);
    double ct3[4] = CROSS(&B[m], gradEdotB);
    double ct2[4] = {0, S[1] * divB, S[2] * divB, S[3] * divB};

    double gamma1 = 0.0;
    double gamma2 = 1.0;
    double gamma3 = 1.0;

    //printf("%f %f %f %f\n", B2, B[m+1], B[m+2], B[m+3]);

    /* Current evaluation */
    ffe_sim_ohms_law(sim,
		     &E[m], rotE, divE, 
		     &B[m], rotB, divB, J);

    for (int d=1; d<=3; ++d) {

      /* Maxwell's equations */
      dtE[m+d] = +rotB[d] - J[d];
      dtB[m+d] = -rotE[d];

      /* Addition of hyperbolicity terms */
      if (0) {
	dtE[m+d] += -gamma1 / B2 * ct1[d];
	dtB[m+d] += -gamma2 / B2 * ct2[d] - gamma3 / B2 * ct3[d];
      }
    }
  }


  /* ===========================================================================
   * Fill in the n-th (n=0,1,2,3) field register
   * ===========================================================================
   */
  FOR_ALL_INTERIOR(Ni, Nj, Nk) {

    int m = IND(i,j,k,0);

    E[m+1] = E0[m+1] + dt * dtE[m+1];
    E[m+2] = E0[m+2] + dt * dtE[m+2];
    E[m+3] = E0[m+3] + dt * dtE[m+3];

    B[m+1] = B0[m+1] + dt * dtB[m+1];
    B[m+2] = B0[m+2] + dt * dtB[m+2];
    B[m+3] = B0[m+3] + dt * dtB[m+3];
  }

  cow_dfield_syncguard(sim->electric[RKstep]);
  cow_dfield_syncguard(sim->magnetic[RKstep]);
  cow_dfield_syncguard(sim->electric[RKstep+4]);
  cow_dfield_syncguard(sim->magnetic[RKstep+4]);
}



/*
 * Average Runge-Kutta substeps to complete a full time step
 * =====================================================================
 */
void ffe_sim_average_rk(struct ffe_sim *sim)
{
  int Ni = cow_domain_getnumlocalzonesinterior(sim->domain, 0);
  int Nj = cow_domain_getnumlocalzonesinterior(sim->domain, 1);
  int Nk = cow_domain_getnumlocalzonesinterior(sim->domain, 2);

  int si = cow_dfield_getstride(sim->electric[0], 0);
  int sj = cow_dfield_getstride(sim->electric[0], 1);
  int sk = cow_dfield_getstride(sim->electric[0], 2);

  double *E = cow_dfield_getdatabuffer(sim->electric[0]);
  double *B = cow_dfield_getdatabuffer(sim->magnetic[0]);

  double *dtE0 = cow_dfield_getdatabuffer(sim->electric[4]);
  double *dtB0 = cow_dfield_getdatabuffer(sim->magnetic[4]);
  double *dtE1 = cow_dfield_getdatabuffer(sim->electric[5]);
  double *dtB1 = cow_dfield_getdatabuffer(sim->magnetic[5]);
  double *dtE2 = cow_dfield_getdatabuffer(sim->electric[6]);
  double *dtB2 = cow_dfield_getdatabuffer(sim->magnetic[6]);
  double *dtE3 = cow_dfield_getdatabuffer(sim->electric[7]);
  double *dtB3 = cow_dfield_getdatabuffer(sim->magnetic[7]);

  double dt = sim->status.time_step;


  /* ===========================================================================
   * Average the RK substeps, write result into register 0
   * ===========================================================================
   */
  FOR_ALL_INTERIOR(Ni, Nj, Nk) {

    int m = IND(i,j,k,0);

    E[m+1] += dt/6 * (dtE0[m+1] + 2*dtE1[m+1] + 2*dtE2[m+1] + dtE3[m+1]);
    E[m+2] += dt/6 * (dtE0[m+2] + 2*dtE1[m+2] + 2*dtE2[m+2] + dtE3[m+2]);
    E[m+3] += dt/6 * (dtE0[m+3] + 2*dtE1[m+3] + 2*dtE2[m+3] + dtE3[m+3]);

    B[m+1] += dt/6 * (dtB0[m+1] + 2*dtB1[m+1] + 2*dtB2[m+1] + dtB3[m+1]);
    B[m+2] += dt/6 * (dtB0[m+2] + 2*dtB1[m+2] + 2*dtB2[m+2] + dtB3[m+2]);
    B[m+3] += dt/6 * (dtB0[m+3] + 2*dtB1[m+3] + 2*dtB2[m+3] + dtB3[m+3]);

  }

  cow_dfield_syncguard(sim->electric[0]);
  cow_dfield_syncguard(sim->magnetic[0]);
}



/*
 * Advance the simulation by one full iteration
 * =====================================================================
 */
void ffe_sim_advance(struct ffe_sim *sim)
{
  double dx = sim->grid_spacing[1];
  double dy = sim->grid_spacing[2];
  double dz = sim->grid_spacing[3];
  double dt = sim->cfl_parameter * MIN3(dx,dy,dz);

  sim->status.time_step = dt;

  ffe_sim_advance_rk(sim, 0);
  ffe_sim_advance_rk(sim, 1);
  ffe_sim_advance_rk(sim, 2);
  ffe_sim_advance_rk(sim, 3);
  ffe_sim_average_rk(sim);

  sim->status.iteration += 1;
  sim->status.time_simulation += dt;
}



int main(int argc, char **argv)
{
  cow_init(0, NULL, 0);

  struct ffe_sim sim;

  ffe_sim_init(&sim);
  ffe_sim_initial_data(&sim);

  int local_grid_size = cow_domain_getnumlocalzonesinterior(sim.domain,
							    COW_ALL_DIMS);

  cow_dfield_write(sim.magnetic[0], "chkpt.0000.h5");
  cow_dfield_write(sim.electric[0], "chkpt.0000.h5");

  while (sim.status.time_simulation < 0.1) {

    clock_t start_cycle, stop_cycle;
    start_cycle = clock();

    ffe_sim_advance(&sim);

    stop_cycle = clock();
    sim.status.kzps = 1e-3 * local_grid_size / (stop_cycle - start_cycle) *
      CLOCKS_PER_SEC;

    printf("n=%06d t=%3.2e %3.2f kzps\n",
	   sim.status.iteration,
	   sim.status.time_simulation,
	   sim.status.kzps);
  }

  cow_dfield_write(sim.magnetic[0], "chkpt.0001.h5");
  cow_dfield_write(sim.electric[0], "chkpt.0001.h5");


  ffe_sim_free(&sim);

  cow_finalize();
  return 0;
}










/*
 * Scratch -------------------------------------------


double GE[4][4] = {0, 0, 0, 0,
		   0, D1(E,1), D1(E,2), D1(E,3),
		   0, D2(E,1), D2(E,2), D2(E,3),
		   0, D3(E,1), D3(E,2), D3(E,3)};
double GB[4][4] = {0, 0, 0, 0,
		   0, D1(B,1), D1(B,2), D1(B,3),
		   0, D2(B,1), D2(B,2), D2(B,3),
		   0, D3(B,1), D3(B,2), D3(B,3)};

*/
