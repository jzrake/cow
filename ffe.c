
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <complex.h>
#include "jsw_rand.h"
#include "cow.h"

#ifndef M_PI
#define M_PI 3.1415926535897926
#endif



#define FFE_NG 3 /* number of guard zones */
#define FFE_DIFFERENCE_ORDER 4
#define FFE_DISSIPATION_ORDER 4

/*
 * Macro for a three-dimensional loop over all interior cells
 * =====================================================================
 */
#define FOR_ALL_INTERIOR(N1, N2, N3)				\
  for (int i=N1==1?0:FFE_NG; i<N1+(N1==1?0:FFE_NG); ++i)	\
    for (int j=N2==1?0:FFE_NG; j<N2+(N2==1?0:FFE_NG); ++j)	\
      for (int k=N3==1?0:FFE_NG; k<N3+(N3==1?0:FFE_NG); ++k)	\

/*
 * Macro to calculate linear index of (i,j,k,m) ... m goes from 1, not 0
 * =====================================================================
 */
#define IND(i,j,k,m) ((i)*si + (j)*sj + (k)*sk) + (m-1)

/*
 * Macros to calculate finite differences
 * =====================================================================
 */
#define DIFF1C2(F,s) ((-1*(F)[-1*s] +		\
		       +1*(F)[+1*s]) / 2.0)

#define DIFF1C4(F,s) ((+1*(F)[-2*s] +		\
		       -8*(F)[-1*s] +		\
		       +8*(F)[+1*s] +		\
		       -1*(F)[+2*s]) / 12.0)

#define DIFF2C2(F,s) ((+1*(F)[-1*s] +		\
		       -2*(F)[+0*s] +		\
		       +1*(F)[+1*s]) / 1.0)

#define DIFF2C4(F,s) ((-1 *(F)[-2*s] +		\
		       +16*(F)[-1*s] +		\
		       -30*(F)[+0*s] +		\
		       +16*(F)[+1*s] +		\
		       -1 *(F)[+2*s]) / 12.0)

#define DIFF4C2(F,s) ((+1*(F)[-2*s] +		\
		       -4*(F)[-1*s] +		\
		       +6*(F)[+0*s] +		\
		       -4*(F)[+1*s] +		\
		       +1*(F)[+2*s]) / 1.0)

#define DIFF4C4(F,s) ((-1 *(F)[-3*s] +		\
		       +12*(F)[-2*s] +		\
		       -39*(F)[-1*s] +		\
		       +56*(F)[+0*s] +		\
		       -39*(F)[+1*s] +		\
		       +12*(F)[+2*s] +		\
		       -1 *(F)[+3*s]) / 6.0)

#define DIFF6C2(F,s) (( 1 *(F)[-3*s] +		\
		       -6 *(F)[-2*s] +		\
		       +15*(F)[-1*s] +		\
		       -20*(F)[+0*s] +		\
		       +15*(F)[+1*s] +		\
		       -6 *(F)[+2*s] +		\
		        1 *(F)[+3*s]) / 1.0)



#define CROSS(E,B) {0.0,				\
      (E)[2]*(B)[3]-(E)[3]*(B)[2],			\
      (E)[3]*(B)[1]-(E)[1]*(B)[3],			\
      (E)[1]*(B)[2]-(E)[2]*(B)[1]}			\

#define DOT(E,B) ((E)[1]*(B)[1] + (E)[2]*(B)[2] + (E)[3]*(B)[3])

#define MAX3(a,b,c) (a>b && b>c ? a : (b > c ? b : c))
#define MIN3(a,b,c) (a<b && b<c ? a : (b < c ? b : c))




struct ffe_sim;

/*
 * Initial data library
 * =====================================================================
 */
typedef void (*InitialDataFunction)(struct ffe_sim *sim, double x[4], double E[4], double B[4]);
static void initial_data_emwave    (struct ffe_sim *sim, double x[4], double E[4], double B[4]);
static void initial_data_alfvenwave(struct ffe_sim *sim, double x[4], double E[4], double B[4]);
static void initial_data_abc       (struct ffe_sim *sim, double x[4], double E[4], double B[4]);
static void initial_data_beltrami  (struct ffe_sim *sim, double x[4], double E[4], double B[4]);
static void random_beltrami_field(double x[4], double B[4], int model, int k2);


enum FfeSimParameter {
  FFE_OHMS_LAW_VACUUM,
  FFE_OHMS_LAW_FORCE_FREE
} ;


struct ffe_measure
{
  double electric_energy;
  double magnetic_energy;
  double magnetic_helicity;
  double magnetic_monopole; /* |div(B)| */
} ;


struct ffe_status
{
  int iteration;
  int checkpoint_number;
  double time_simulation;
  double time_final;
  double time_step;
  double time_last_checkpoint;
  double kzps;
} ;


struct ffe_sim
{
  cow_domain *domain;
  cow_dfield *electric[6]; /* 0,1: E, 4-6: dtE */
  cow_dfield *magnetic[6];

  int Ni, Nj, Nk;
  double grid_spacing[4];
  double cfl_parameter; /* Courant number [0.0 - 0.25]*/
  double eps_parameter; /* Kreiss-Oliger parameter [0 - 1] */

  enum FfeSimParameter flag_ohms_law;
  InitialDataFunction initial_data;
  double time_between_checkpoints;
  double time_final;

  struct ffe_status status;

  char output_directory[1024];
} ;



void ffe_sim_init(struct ffe_sim *sim);
void ffe_sim_free(struct ffe_sim *sim);
void ffe_sim_analyze(struct ffe_sim *sim, char *filename);
void ffe_sim_initial_data(struct ffe_sim *sim);
int  ffe_sim_problem_setup(struct ffe_sim *sim, const char *problem_name);



/*
 * Initialze a new simulation instance. Parameters that need to be initialized
 * before this call:
 *
 * -> Ni, Nj, Nk
 * -> initial_data
 * -> flag_ohms_law
 * -> time_between_checkpoints
 * -> time_final
 * -> output_directory
 *
 * =====================================================================
 */
void ffe_sim_init(struct ffe_sim *sim)
{
  sim->domain = cow_domain_new();


  sim->grid_spacing[0] = 0.0;
  sim->grid_spacing[1] = 1.0 / sim->Ni;
  sim->grid_spacing[2] = 1.0 / sim->Nj;
  sim->grid_spacing[3] = 1.0 / sim->Nk;
  sim->cfl_parameter = 0.25;
  sim->eps_parameter = 1.00; /* [0-1] */

  sim->status.iteration = 0;
  sim->status.checkpoint_number = 0;
  sim->status.time_simulation = 0.0;
  sim->status.time_step = 0.0;
  sim->status.time_last_checkpoint = -sim->time_between_checkpoints;



  cow_domain_setndim(sim->domain, (sim->Ni>1) + (sim->Nj>1) + (sim->Nk>1));
  cow_domain_setsize(sim->domain, 0, sim->Ni);
  cow_domain_setsize(sim->domain, 1, sim->Nj);
  cow_domain_setsize(sim->domain, 2, sim->Nk);
  cow_domain_setguard(sim->domain, FFE_NG);
  cow_domain_commit(sim->domain);

  for (int n=0; n<6; ++n) {
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

  for (int n=0; n<6; ++n) {
    cow_dfield_del(sim->electric[n]);
    cow_dfield_del(sim->magnetic[n]);
  }

  cow_domain_del(sim->domain);
}



int ffe_sim_problem_setup(struct ffe_sim *sim, const char *problem_name)
{
  if (problem_name == NULL) {
    printf("1. emwave\n");
    printf("2. alfvenwave\n");
    printf("3. abc\n");
    printf("4. beltrami\n");
    return 0;
  }
  else if (!strcmp(problem_name, "emwave")) {
    sim->Ni = 16;
    sim->Nj = 16;
    sim->Nk = 64;
    sim->initial_data = initial_data_emwave;
    sim->flag_ohms_law = FFE_OHMS_LAW_VACUUM;
    return 0;
  }
  else if (!strcmp(problem_name, "alfvenwave")) {
    sim->Ni = 16;
    sim->Nj = 16;
    sim->Nk = 64;
    sim->initial_data = initial_data_alfvenwave;
    sim->flag_ohms_law = FFE_OHMS_LAW_FORCE_FREE;
    return 0;
  }
  else if (!strcmp(problem_name, "abc")) {
    sim->Ni = 128;
    sim->Nj = 128;
    sim->Nk = 1;
    sim->initial_data = initial_data_abc;
    sim->flag_ohms_law = FFE_OHMS_LAW_FORCE_FREE;
    return 0;
  }
  else if (!strcmp(problem_name, "beltrami")) {
    sim->Ni = 64;
    sim->Nj = 64;
    sim->Nk = 64;
    sim->initial_data = initial_data_beltrami;
    sim->flag_ohms_law = FFE_OHMS_LAW_FORCE_FREE;
    return 0;
  }
  else {
    return 1;
  }
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
  int si = cow_dfield_getstride(sim->electric[0], 0);
  int sj = cow_dfield_getstride(sim->electric[0], 1);
  int sk = cow_dfield_getstride(sim->electric[0], 2);
  double *E = cow_dfield_getdatabuffer(sim->electric[0]);
  double *B = cow_dfield_getdatabuffer(sim->magnetic[0]);

  FOR_ALL_INTERIOR(Ni, Nj, Nk) {    

    int m = IND(i,j,k,0);
    double x[4] = {0,
    		   cow_domain_positionatindex(sim->domain, 0, i),
    		   cow_domain_positionatindex(sim->domain, 1, j),
    		   cow_domain_positionatindex(sim->domain, 2, k)};

    sim->initial_data(sim, x, &E[m], &B[m]);
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

      if (B2 < 1e-12) B2 = 1e-12; /* Don't divide by zero! */

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
/* https://einsteintoolkit.org/documentation/ThornDoc/CactusNumerical/Dissipation */
#if (FFE_DIFFERENCE_ORDER == 2)
#define D1(F,c)  (Ni==1 ? 0.0 : DIFF1C2(F+m+c,si)/dx)
#define D2(F,c)  (Nj==1 ? 0.0 : DIFF1C2(F+m+c,sj)/dy)
#define D3(F,c)  (Nk==1 ? 0.0 : DIFF1C2(F+m+c,sk)/dz)
#elif (FFE_DIFFERENCE_ORDER == 4)
#define D1(F,c)  (Ni==1 ? 0.0 : DIFF1C4(F+m+c,si)/dx)
#define D2(F,c)  (Nj==1 ? 0.0 : DIFF1C4(F+m+c,sj)/dy)
#define D3(F,c)  (Nk==1 ? 0.0 : DIFF1C4(F+m+c,sk)/dz)
#endif
#if (FFE_DISSIPATION_ORDER == 2)
#define KO(F,c) ((Ni==1 ? 0.0 : DIFF2C2(F+m+c,si)/dt) +    \
                 (Nj==1 ? 0.0 : DIFF2C2(F+m+c,sj)/dt) +    \
                 (Nk==1 ? 0.0 : DIFF2C2(F+m+c,sk)/dt))
  double KO_const = -1.0/4;
#elif (FFE_DISSIPATION_ORDER == 4)
#define KO(F,c) ((Ni==1 ? 0.0 : DIFF4C2(F+m+c,si)/dt) +	   \
		 (Nj==1 ? 0.0 : DIFF4C2(F+m+c,sj)/dt) +	   \
                 (Nk==1 ? 0.0 : DIFF4C2(F+m+c,sk)/dt))
  double KO_const = +1.0/16;
#elif (FFE_DISSIPATION_ORDER == 6)
#define KO(F,c) ((Ni==1 ? 0.0 : DIFF6C2(F+m+c,si)/dt) +	   \
		 (Nj==1 ? 0.0 : DIFF6C2(F+m+c,sj)/dt) +	   \
                 (Nk==1 ? 0.0 : DIFF6C2(F+m+c,sk)/dt))
  double KO_const = -1.0/64;
#endif


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
  double *E  = cow_dfield_getdatabuffer(sim->electric[1]);
  double *B  = cow_dfield_getdatabuffer(sim->magnetic[1]);

  double *dtE = cow_dfield_getdatabuffer(sim->electric[RKstep+2]);
  double *dtB = cow_dfield_getdatabuffer(sim->magnetic[RKstep+2]);

  double dx = sim->grid_spacing[1];
  double dy = sim->grid_spacing[2];
  double dz = sim->grid_spacing[3];
  double dt = sim->status.time_step;

  /* ===========================================================================
   * Fill in the n-th (n=0,1,2,3) time-derivative register, reading from the
   * [1] field register
   * ===========================================================================
   */
  FOR_ALL_INTERIOR(Ni, Nj, Nk) {


    int m = IND(i,j,k,0);
    double divE = D1(E,1) + D2(E,2) + D3(E,3);
    double divB = D1(B,1) + D2(B,2) + D3(B,3);    
    double rotE[4] = {0, D2(E,3) - D3(E,2), D3(E,1) - D1(E,3), D1(E,2) - D2(E,1)};
    double rotB[4] = {0, D2(B,3) - D3(B,2), D3(B,1) - D1(B,3), D1(B,2) - D2(B,1)};


    /* Current evaluation */
    double J[4];
    ffe_sim_ohms_law(sim, &E[m], rotE, divE, &B[m], rotB, divB, J);


    /* Maxwell's equations */
    for (int d=1; d<=3; ++d) {
      dtE[m+d] = +rotB[d] - J[d];
      dtB[m+d] = -rotE[d];
    }



    /* Kreiss-Oliger dissipation */
    if (0) {
      double lplE[4] = {0, KO(E,1), KO(E,2), KO(E,3)};
      double lplB[4] = {0, KO(B,1), KO(B,2), KO(B,3)};
      double eps = sim->eps_parameter;

      for (int d=1; d<=3; ++d) {
	dtE[m+d] -= eps * KO_const * lplE[d];
	dtB[m+d] -= eps * KO_const * lplB[d];
      }
    }



    /* Hyperbolicity terms, eqn 48-49: Pfeiffer (2013) */
    if (0) {
      double B2 = DOT(&B[m], &B[m]);
      double gradEdotB[4] = {0, 0, 0, 0};

      double S[4]   = CROSS(&E[m], &B[m]);
      double ct1[4] = CROSS(&E[m], gradEdotB);
      double ct3[4] = CROSS(&B[m], gradEdotB);
      double ct2[4] = {0, S[1] * divB, S[2] * divB, S[3] * divB};

      double gamma1 = 0.0;
      double gamma2 = 1.0;
      double gamma3 = 1.0;

      if (B2 < 1e-12) B2 = 1e-12;

      for (int d=1; d<=3; ++d) {
	dtE[m+d] -= gamma1 / B2 * ct1[d];
	dtB[m+d] -= gamma2 / B2 * ct2[d] + gamma3 / B2 * ct3[d];
      }
    }
  }



  /* ===========================================================================
   * Fill in the n-th (n=0,1,2,3) field register
   * ===========================================================================
   */
  FOR_ALL_INTERIOR(Ni, Nj, Nk) {

    int m = IND(i,j,k,0);

    E[m+1] = E0[m+1] + dt * RKparam * dtE[m+1];
    E[m+2] = E0[m+2] + dt * RKparam * dtE[m+2];
    E[m+3] = E0[m+3] + dt * RKparam * dtE[m+3];

    B[m+1] = B0[m+1] + dt * RKparam * dtB[m+1];
    B[m+2] = B0[m+2] + dt * RKparam * dtB[m+2];
    B[m+3] = B0[m+3] + dt * RKparam * dtB[m+3];

  }

  cow_dfield_syncguard(sim->electric[0]);
  cow_dfield_syncguard(sim->magnetic[0]);
  cow_dfield_syncguard(sim->electric[1]);
  cow_dfield_syncguard(sim->magnetic[1]);

#undef D1
#undef D2
#undef D3
#undef KO
}



/*
 * Apply Kreiss-Oliger operator to subtract high frequencies
 * =====================================================================
 */
void ffe_sim_kreiss_oliger(struct ffe_sim *sim)
{
#if (FFE_DISSIPATION_ORDER == 2)
#define KO(F,c) ((Ni==1 ? 0.0 : DIFF2C2(F+m+c,si)) +	\
		 (Nj==1 ? 0.0 : DIFF2C2(F+m+c,sj)) +	\
		 (Nk==1 ? 0.0 : DIFF2C2(F+m+c,sk)))
  double KO_const = -1.0/4;
#elif (FFE_DISSIPATION_ORDER == 4)
#define KO(F,c) ((Ni==1 ? 0.0 : DIFF4C2(F+m+c,si)) +	\
		 (Nj==1 ? 0.0 : DIFF4C2(F+m+c,sj)) +	\
		 (Nk==1 ? 0.0 : DIFF4C2(F+m+c,sk)))
  double KO_const = +1.0/16;
#elif (FFE_DISSIPATION_ORDER == 6)
#define KO(F,c) ((Ni==1 ? 0.0 : DIFF6C2(F+m+c,si)) +	   \
		 (Nj==1 ? 0.0 : DIFF6C2(F+m+c,sj)) +	   \
                 (Nk==1 ? 0.0 : DIFF6C2(F+m+c,sk)))
  double KO_const = -1.0/64;
#endif

  int Ni = cow_domain_getnumlocalzonesinterior(sim->domain, 0);
  int Nj = cow_domain_getnumlocalzonesinterior(sim->domain, 1);
  int Nk = cow_domain_getnumlocalzonesinterior(sim->domain, 2);
  int si = cow_dfield_getstride(sim->electric[0], 0);
  int sj = cow_dfield_getstride(sim->electric[0], 1);
  int sk = cow_dfield_getstride(sim->electric[0], 2);
  double *E = cow_dfield_getdatabuffer(sim->electric[0]);
  double *B = cow_dfield_getdatabuffer(sim->magnetic[0]);
  double *dE = cow_dfield_getdatabuffer(sim->electric[2]); /* use dt registers */
  double *dB = cow_dfield_getdatabuffer(sim->magnetic[2]);


  /* ===========================================================================
   * Fill in the n-th (n=0,1,2,3) time-derivative register, reading from the
   * [1] field register
   * ===========================================================================
   */
  FOR_ALL_INTERIOR(Ni, Nj, Nk) {

    int m = IND(i,j,k,0);

    for (int d=1; d<=3; ++d) {
      dE[m+d] = KO(E,d);
      dB[m+d] = KO(B,d);
    }

  }

  FOR_ALL_INTERIOR(Ni, Nj, Nk) {

    int m = IND(i,j,k,0);
    double eps = sim->eps_parameter;

    for (int d=1; d<=3; ++d) {
      E[m+d] -= eps * KO_const * dE[m+d];
      B[m+d] -= eps * KO_const * dB[m+d];
    }

  }


  cow_dfield_syncguard(sim->electric[0]);
  cow_dfield_syncguard(sim->magnetic[0]);

#undef KO
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

  double *dtE0 = cow_dfield_getdatabuffer(sim->electric[2]);
  double *dtB0 = cow_dfield_getdatabuffer(sim->magnetic[2]);
  double *dtE1 = cow_dfield_getdatabuffer(sim->electric[3]);
  double *dtB1 = cow_dfield_getdatabuffer(sim->magnetic[3]);
  double *dtE2 = cow_dfield_getdatabuffer(sim->electric[4]);
  double *dtB2 = cow_dfield_getdatabuffer(sim->magnetic[4]);
  double *dtE3 = cow_dfield_getdatabuffer(sim->electric[5]);
  double *dtB3 = cow_dfield_getdatabuffer(sim->magnetic[5]);

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




    /*
     * Subtract out the component of E parallel to B
     */
    double BB = DOT(&B[m], &B[m]);
    double EB = DOT(&E[m], &B[m]);


    E[m+1] -= EB/BB * B[m+1];
    E[m+2] -= EB/BB * B[m+2];
    E[m+3] -= EB/BB * B[m+3];




    /*
     * Cap the electric field vector to ensure E <= B
     */
    double EE = DOT(&E[m], &E[m]);

    if (EE > BB) {

      double f = sqrt(BB/EE);

      E[m+1] *= f;
      E[m+2] *= f;
      E[m+3] *= f;

    }

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

  /* copy data from [0] register into [1] register */
  cow_dfield_copy(sim->electric[0], sim->electric[1]);
  cow_dfield_copy(sim->magnetic[0], sim->magnetic[1]);

  ffe_sim_advance_rk(sim, 0);
  ffe_sim_advance_rk(sim, 1);
  ffe_sim_advance_rk(sim, 2);
  ffe_sim_advance_rk(sim, 3);
  ffe_sim_average_rk(sim);
  ffe_sim_kreiss_oliger(sim);

  sim->status.iteration += 1;
  sim->status.time_simulation += dt;
}



/*
 * Carry out measurement diagnostic
 * =====================================================================
 */
void ffe_sim_measure(struct ffe_sim *sim, struct ffe_measure *meas)
{
#if (FFE_DIFFERENCE_ORDER == 2)

#define D1(F,c)  (Ni==1 ? 0.0 : DIFF1C2(F+m+c,si)/dx)
#define D2(F,c)  (Nj==1 ? 0.0 : DIFF1C2(F+m+c,sj)/dy)
#define D3(F,c)  (Nk==1 ? 0.0 : DIFF1C2(F+m+c,sk)/dz)

#elif (FFE_DIFFERENCE_ORDER == 4)

#define D1(F,c)  (Ni==1 ? 0.0 : DIFF1C4(F+m+c,si)/dx)
#define D2(F,c)  (Nj==1 ? 0.0 : DIFF1C4(F+m+c,sj)/dy)
#define D3(F,c)  (Nk==1 ? 0.0 : DIFF1C4(F+m+c,sk)/dz)

#endif

#define GLB_AVG(x) x = cow_domain_dblsum(sim->domain, x) / Nt

  int Ni = cow_domain_getnumlocalzonesinterior(sim->domain, 0);
  int Nj = cow_domain_getnumlocalzonesinterior(sim->domain, 1);
  int Nk = cow_domain_getnumlocalzonesinterior(sim->domain, 2);
  int Nt = cow_domain_getnumglobalzones(sim->domain, COW_ALL_DIMS);
  int si = cow_dfield_getstride(sim->electric[0], 0);
  int sj = cow_dfield_getstride(sim->electric[0], 1);
  int sk = cow_dfield_getstride(sim->electric[0], 2);
  double dx = sim->grid_spacing[1];
  double dy = sim->grid_spacing[2];
  double dz = sim->grid_spacing[3];

  double *E = cow_dfield_getdatabuffer(sim->electric[0]);
  double *B = cow_dfield_getdatabuffer(sim->magnetic[0]);

  meas->electric_energy = 0.0;
  meas->magnetic_energy = 0.0;
  meas->magnetic_monopole = 0.0; /* TODO */
  meas->magnetic_helicity = 0.0;

  FOR_ALL_INTERIOR(Ni, Nj, Nk) {

    int m = IND(i,j,k,0);

    double divB = D1(B,1) + D2(B,2) + D3(B,3);    

    double EE = DOT(&E[m], &E[m]);
    double BB = DOT(&B[m], &B[m]);

    meas->electric_energy += 0.5 * EE;
    meas->magnetic_energy += 0.5 * BB;
    meas->magnetic_monopole += fabs(divB);

  }

  GLB_AVG(meas->electric_energy);
  GLB_AVG(meas->magnetic_energy);
  GLB_AVG(meas->magnetic_monopole);

#undef D1
#undef D2
#undef D3
}



int main(int argc, char **argv)
{
  cow_init(0, NULL, 0);



  const char *problem_name = NULL;
  char logfile_name[1024];
  struct ffe_sim sim;
  struct ffe_measure measure;

  sim.time_final = 1.0;
  sim.time_between_checkpoints = 1.0;
  strcpy(sim.output_directory, ".");



  /*
   * Print a help message
   * ===================================================================
   */

  printf("\nForce-free electrodynamics solver\n");
  printf("Jonathan Zrake, Stanford University (2015)\n");

  if (argc == 1) {
    printf("usage: ffe <problem-name> [sim.time_final=1.0] [N=16,16,16]\n");
    printf("problems are:\n");
    ffe_sim_problem_setup(NULL, NULL);
    return 0;
  }
  else {
    problem_name = argv[1];
  }



  /*
   * Set up the problem defaults
   * ===================================================================
   */
  if (ffe_sim_problem_setup(&sim, problem_name)) {
    printf("[ffe] error: unkown problem name: '%s', choose one of\n", problem_name);
    ffe_sim_problem_setup(NULL, NULL);
    return 1;
  }



  /*
   * Scan command line arguments
   * ===================================================================
   */
  for (int n=2; n<argc; ++n) {

    if (!strncmp(argv[n], "tmax=", 5)) {
      sscanf(argv[n], "tmax=%lf", &sim.time_final);
    }
    else if (!strncmp(argv[n], "cpi=", 4)) {
      sscanf(argv[n], "cpi=%lf", &sim.time_between_checkpoints);
    }
    else if (!strncmp(argv[n], "N=", 2)) {
      int num = sscanf(argv[n], "N=%d,%d,%d", &sim.Ni, &sim.Nj, &sim.Nk);
      if (num != 3) {
	printf("[ffe] error: badly formatted option '%s'\n", argv[n]);
	return 1;
      }
    }
    else if (!strncmp(argv[n], "outdir=", 7)) {
      sscanf(argv[n], "outdir=%1024s", sim.output_directory);
    }
    else if (!strncmp(argv[n], "eps=", 4)) {
      sscanf(argv[n], "eps=%lf", &sim.eps_parameter);
    }
    else {
      printf("[ffe] error: unrecognized option '%s'\n", argv[n]);
      return 1;
    }
  }





  ffe_sim_init(&sim);
  ffe_sim_initial_data(&sim);


  printf("\n-----------------------------------------\n");
  printf("tmax ....................... %12.10lf\n", sim.time_final);
  printf("time_between_checkpoints ... %12.10lf\n", sim.time_between_checkpoints);
  printf("eps_parameter .............. %12.10lf\n", sim.eps_parameter);
  printf("output_directory ........... %s\n", sim.output_directory);
  printf("-----------------------------------------\n\n");


  int invalid_output = 0;
  int local_grid_size = cow_domain_getnumlocalzonesinterior(sim.domain,
							    COW_ALL_DIMS);


  snprintf(logfile_name, 1024, "%s/ffe.dat", sim.output_directory);

  if (cow_domain_getcartrank(sim.domain) == 0) {
  
    FILE *logf = fopen(logfile_name, "w");
    
    if (logf == NULL) {
      printf("[ffe] error: could not open log file '%s'\n", logfile_name);
      invalid_output = 1;
    }
    else {
      fclose(logf);
      invalid_output = 0;
    }
  }


  /* Propagate any errors to all procs */
  invalid_output = cow_domain_intsum(sim.domain, invalid_output);  

  if (invalid_output) {
    sim.time_final = 0.0;
  }


  while (sim.status.time_simulation < sim.time_final) {


    /*
     * Write a checkpoint if it's time
     * =================================================================
     */
    if (sim.status.time_simulation - sim.status.time_last_checkpoint >=
	sim.time_between_checkpoints) {

      char chkpt_name[1024];
      snprintf(chkpt_name, 1024, "%s/chkpt.%04d.h5",
	       sim.output_directory,
	       sim.status.checkpoint_number);

      cow_dfield_write(sim.magnetic[0], chkpt_name);
      cow_dfield_write(sim.electric[0], chkpt_name);

      sim.status.time_last_checkpoint += sim.time_between_checkpoints;
      sim.status.checkpoint_number += 1;
    }


    /*
     * Handle post-processing and reductions
     * =================================================================
     */

    ffe_sim_measure(&sim, &measure);

    if (sim.status.iteration % 100 == 0) {
      ffe_sim_analyze(&sim, "analysis.h5");
    }

    if (cow_domain_getcartrank(sim.domain) == 0) {

      FILE *logf = fopen(logfile_name, "a");

      fprintf(logf, "%+12.10e %+12.10e %+12.10e %+12.10e %+12.10e\n",
	      sim.status.time_simulation,
	      measure.electric_energy,
	      measure.magnetic_energy,
	      measure.magnetic_helicity,
	      measure.magnetic_monopole);

      fclose(logf);
    }



    /*
     * Evolve the system
     * =================================================================
     */
    clock_t start_cycle, stop_cycle;
    start_cycle = clock();

    ffe_sim_advance(&sim);

    stop_cycle = clock();
    sim.status.kzps = 1e-3 * local_grid_size / (stop_cycle - start_cycle) *
      CLOCKS_PER_SEC;

    printf("[ffe] n=%06d t=%3.2e %3.2f kzps\n",
	   sim.status.iteration,
	   sim.status.time_simulation,
	   sim.status.kzps);
  }



  if (invalid_output == 0) {
    char chkpt_name[1024];
    snprintf(chkpt_name, 1024, "%s/chkpt.%04d.h5",
	     sim.output_directory,
	     sim.status.checkpoint_number);
    cow_dfield_write(sim.magnetic[0], chkpt_name);
    cow_dfield_write(sim.electric[0], chkpt_name);
  }


  ffe_sim_free(&sim);

  cow_finalize();
  return 0;
}








void initial_data_emwave(struct ffe_sim *sim, double x[4], double E[4], double B[4])
{
  E[1] = sin(2 * M_PI * x[3]);
  E[2] = 0.0;
  E[3] = 0.0;

  B[1] = 0.0;
  B[2] = sin(2 * M_PI * x[3]);
  B[3] = 0.0;
}

void initial_data_alfvenwave(struct ffe_sim *sim, double x[4], double E[4], double B[4])
{
  E[1] = 0.0;
  E[2] = 0.0;
  E[3] = 0.0;

  B[1] = 0.0;
  B[2] = 0.1 * sin(2 * M_PI * x[3]);
  B[3] = 1.0;
}

void initial_data_abc(struct ffe_sim *sim, double x[4], double E[4], double B[4])
{
  double a=1, b=1, c=0;
  double alpha = 2 * M_PI * 2;

  E[1] = 0.0;
  E[2] = 0.0;
  E[3] = 0.0;

  B[1] = 0.0;
  B[2] = 0.0;
  B[3] = 0.0;

  B[2] += a * cos(alpha * x[1]);
  B[3] -= a * sin(alpha * x[1]);
  B[1] += 0.0;

  B[3] += b * cos(alpha * x[2]);
  B[1] -= b * sin(alpha * x[2]);
  B[2] += 0.0;

  B[1] += c * cos(alpha * x[3]);
  B[2] -= c * sin(alpha * x[3]);
  B[3] += 0.0;
}

void initial_data_beltrami(struct ffe_sim *sim, double x[4], double E[4], double B[4])
{
  E[1] = 0.0;
  E[2] = 0.0;
  E[3] = 0.0;

  double X[4] = {0,
		 -1 + 2*x[1],
		 -1 + 2*x[2],
		 -1 + 2*x[3]};

  random_beltrami_field(X, B, 0, 11);
}












































void ffe_sim_analyze(struct ffe_sim *sim, char *filename)
{
  cow_domain *domain = sim->domain;
  cow_dfield *magnetic = sim->magnetic[0];
  cow_dfield *electric = sim->electric[0];
  cow_dfield *vecpoten = cow_dfield_new();
  cow_dfield *jcurrent = cow_dfield_new();



  /* Data fields setup */
  /* ---------------------------------------------------- */
  cow_dfield_setdomain(vecpoten, domain);
  cow_dfield_setname(vecpoten, "vector_potential");
  cow_dfield_addmember(vecpoten, "A1");
  cow_dfield_addmember(vecpoten, "A2");
  cow_dfield_addmember(vecpoten, "A3");
  cow_dfield_commit(vecpoten);

  cow_dfield_setdomain(jcurrent, domain);
  cow_dfield_setname(jcurrent, "electric_current");
  cow_dfield_addmember(jcurrent, "J1");
  cow_dfield_addmember(jcurrent, "J2");
  cow_dfield_addmember(jcurrent, "J3");
  cow_dfield_commit(jcurrent);



  /* Histograms setup */
  /* ---------------------------------------------------- */
  cow_histogram *Pb = cow_histogram_new();
  cow_histogram_setlower(Pb, 0, 1);
  cow_histogram_setupper(Pb, 0, 8192);
  cow_histogram_setnbins(Pb, 0, 4096);
  cow_histogram_setspacing(Pb, COW_HIST_SPACING_LINEAR);
  cow_histogram_setfullname(Pb, "magnetic");
  cow_histogram_setnickname(Pb, "magnetic");


  cow_histogram *Pe = cow_histogram_new();
  cow_histogram_setlower(Pe, 0, 1);
  cow_histogram_setupper(Pe, 0, 8192);
  cow_histogram_setnbins(Pe, 0, 4096);
  cow_histogram_setspacing(Pe, COW_HIST_SPACING_LINEAR);
  cow_histogram_setfullname(Pe, "electric");
  cow_histogram_setnickname(Pe, "electric");


  cow_histogram *Hr = cow_histogram_new();
  cow_histogram_setlower(Hr, 0, 1);
  cow_histogram_setupper(Hr, 0, 8192);
  cow_histogram_setnbins(Hr, 0, 4096);
  cow_histogram_setspacing(Hr, COW_HIST_SPACING_LINEAR);
  cow_histogram_setfullname(Hr, "helicity-real");
  cow_histogram_setnickname(Hr, "helicity-real");


  cow_histogram *Hi = cow_histogram_new();
  cow_histogram_setlower(Hi, 0, 1);
  cow_histogram_setupper(Hi, 0, 8192);
  cow_histogram_setnbins(Hi, 0, 4096);
  cow_histogram_setspacing(Hi, COW_HIST_SPACING_LINEAR);
  cow_histogram_setfullname(Hi, "helicity-imag");
  cow_histogram_setnickname(Hi, "helicity-imag");


  cow_fft_inversecurl(magnetic, vecpoten);
  cow_fft_curl       (magnetic, jcurrent);



  cow_fft_pspecvecfield(magnetic, Pb);
  cow_fft_pspecvecfield(electric, Pe);
  cow_fft_helicityspec(magnetic, Hr, Hi);

  if (filename) {

    char gname[1024];

    snprintf(gname, 1024, "spectra-%06d", sim->status.iteration);

    cow_histogram_dumphdf5(Pb, filename, gname);
    cow_histogram_dumphdf5(Pe, filename, gname);
    cow_histogram_dumphdf5(Hi, filename, gname);
    cow_histogram_dumphdf5(Hr, filename, gname);


    if (0) { /* write derived fields */
      cow_dfield_write(magnetic, filename);
      cow_dfield_write(vecpoten, filename);
      cow_dfield_write(jcurrent, filename);
    }

  }

  cow_histogram_del(Pb);
  cow_histogram_del(Pe);
  cow_histogram_del(Hr);
  cow_histogram_del(Hi);

  cow_dfield_del(vecpoten);
  cow_dfield_del(jcurrent);
}









typedef complex double Complex;
typedef struct fourier_mode
{
  double k[4];
  Complex A[4];
} fourier_mode;




void random_beltrami_field(double x[4], double B[4], int model, int k2)
{
#define RAND jsw_random_double(&rand, -1, 1)
  int m,i,j,k;
  jsw_rand_t rand;
  jsw_seed(&rand, model);


  int k2_sphere = k2;
  int k_cube = floor(sqrt(k2_sphere)) + 1;
  fourier_mode *modes = NULL;
  int num_modes = 0;
  Complex A[4] = {0, 0, 0, 0};
  double amp;

  for (i=-k_cube; i<=k_cube; ++i) {
    for (j=-k_cube; j<=k_cube; ++j) {
      for (k=-k_cube; k<=k_cube; ++k) {
  	fourier_mode M;

  	if (i*i + j*j + k*k != k2_sphere) {
  	  continue;
  	}
  	else {
  	  M.k[0] = 0.0;
  	  M.k[1] = i;
  	  M.k[2] = j;
  	  M.k[3] = k;
  	  M.A[0] = 0.0;
  	  M.A[1] = RAND + RAND*I;
  	  M.A[2] = RAND + RAND*I;
  	  M.A[3] = RAND + RAND*I;
  	  amp = sqrt(M.A[1]*conj(M.A[1]) + M.A[2]*conj(M.A[2]) + M.A[3]*conj(M.A[3]));
  	  M.A[1] /= amp;
  	  M.A[2] /= amp;
  	  M.A[3] /= amp;
  	  num_modes += 1;
  	  modes = (fourier_mode *) realloc(modes, num_modes * sizeof(fourier_mode));
  	  modes[num_modes-1] = M;

  	  /* printf("k[%d] = [%d %d %d] is on shell\n", num_modes, i, j, k); */
  	}
      }
    }
  }


  for (m=0; m<num_modes; ++m) {
    fourier_mode M = modes[m];
    double a = sqrt(M.k[1]*M.k[1] + M.k[2]*M.k[2] + M.k[3]*M.k[3]);
    Complex K[4] = {0, I*M.k[1], I*M.k[2], I*M.k[3]};
    Complex Ikx  = (K[1]*x[1] + K[2]*x[2] + K[3]*x[3]) * M_PI;
    Complex P[4] = {0, M.A[1], M.A[2], M.A[3]}; /* a times psi */

    Complex T[4] = {0, /* T = K cross (a psi) */
		    K[2]*P[3] - K[3]*P[2],
		    K[3]*P[1] - K[1]*P[3],
		    K[1]*P[2] - K[2]*P[1]};

    Complex S[4] = {0, /* S = K cross T / alpha */
		    (K[2]*T[3] - K[3]*T[2])/a,
		    (K[3]*T[1] - K[1]*T[3])/a,
		    (K[1]*T[2] - K[2]*T[1])/a};

    A[1] += (S[1] + T[1])/a * cexp(Ikx);
    A[2] += (S[2] + T[2])/a * cexp(Ikx);
    A[3] += (S[3] + T[3])/a * cexp(Ikx);
  }


  free(modes);

  B[1] = A[1];
  B[2] = A[2];
  B[3] = A[3];

#undef RAND
}
