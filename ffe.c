
#include <math.h>
#include <stdio.h>
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
#define GRAD1(F,m) ((F[IND(i+1,j,k,m)] - F[IND(i-1,j,k,m)]) / (2*dx))
#define GRAD2(F,m) ((F[IND(i,j+1,k,m)] - F[IND(i,j-1,k,m)]) / (2*dy))
#define GRAD3(F,m) ((F[IND(i,j,k+1,m)] - F[IND(i,j,k-1,m)]) / (2*dz))


#define MAX3(a,b,c) (a>b && b>c ? a : (b > c ? b : c))
#define MIN3(a,b,c) (a<b && b<c ? a : (b < c ? b : c))

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
} ;


struct ffe_sim
{
  cow_domain *domain;
  cow_dfield *electric[4];
  cow_dfield *magnetic[4];

  int resolution[4];
  double grid_spacing[4];
  double cfl_parameter;

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

  int Ni = 16;
  int Nj = 16;
  int Nk = 64;

  sim->resolution[0] = 0;
  sim->resolution[1] = Ni;
  sim->resolution[2] = Nj;
  sim->resolution[3] = Nk;
  sim->cfl_parameter = 0.1;

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

  for (int n=0; n<4; ++n) {
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
 * Initialze a new simulation instance
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

    E[IND(i,j,k,1)] = sin(2 * M_PI * z);
    E[IND(i,j,k,2)] = 0.0;
    E[IND(i,j,k,3)] = 0.0;

    B[IND(i,j,k,1)] = 0.0;
    B[IND(i,j,k,2)] = cos(2 * M_PI * z);
    B[IND(i,j,k,3)] = 0.0;

  }
}


/*
 * Finalize a simulation instance
 * =====================================================================
 */
void ffe_sim_free(struct ffe_sim *sim)
{

  for (int n=0; n<4; ++n) {
    cow_dfield_del(sim->electric[n]);
    cow_dfield_del(sim->magnetic[n]);
  }

  cow_domain_del(sim->domain);
}


/*
 * Advance the simulation by one iteration
 * =====================================================================
 */
void ffe_sim_advance(struct ffe_sim *sim)
{
  cow_dfield_syncguard(sim->electric[0]);
  cow_dfield_syncguard(sim->magnetic[0]);

  int Ni = cow_domain_getnumlocalzonesinterior(sim->domain, 0);
  int Nj = cow_domain_getnumlocalzonesinterior(sim->domain, 1);
  int Nk = cow_domain_getnumlocalzonesinterior(sim->domain, 2);

  int si = cow_dfield_getstride(sim->electric[0], 0);
  int sj = cow_dfield_getstride(sim->electric[0], 1);
  int sk = cow_dfield_getstride(sim->electric[0], 2);

  double *E = cow_dfield_getdatabuffer(sim->electric[0]);
  double *B = cow_dfield_getdatabuffer(sim->magnetic[0]);

  double *dtE = cow_dfield_getdatabuffer(sim->electric[1]);
  double *dtB = cow_dfield_getdatabuffer(sim->magnetic[1]);

  double dx = sim->grid_spacing[1];
  double dy = sim->grid_spacing[2];
  double dz = sim->grid_spacing[3];

  double dt = sim->cfl_parameter * MIN3(dx,dy,dz);

  //FOR_ALL_INTERIOR(Ni, Nj, Nk) {

  for (int i=2; i<Ni+2; ++i) {
    for (int j=2; j<Nj+2; ++j) {
      for (int k=2; k<Nk+2; ++k) {

	/* double divE = GRAD1(E,1) + GRAD2(E,2) + GRAD3(E,3); */
	/* double divB = GRAD1(B,1) + GRAD2(B,2) + GRAD3(B,3); */

	double rotE1 = GRAD2(E,3) - GRAD3(E,2);
	double rotE2 = GRAD3(E,1) - GRAD1(E,3);
	double rotE3 = GRAD1(E,2) - GRAD2(E,1);

	double rotB1 = GRAD2(B,3) - GRAD3(B,2);
	double rotB2 = GRAD3(B,1) - GRAD1(B,3);
	double rotB3 = GRAD1(B,2) - GRAD2(B,1);

	dtE[IND(i,j,k,1)] = +rotB1;
	dtE[IND(i,j,k,2)] = +rotB2;
	dtE[IND(i,j,k,3)] = +rotB3;

	dtB[IND(i,j,k,1)] = -rotE1;
	dtB[IND(i,j,k,2)] = -rotE2;
	dtB[IND(i,j,k,3)] = -rotE3;
      }
    }
  }

  for (int i=2; i<Ni+2; ++i) {
    for (int j=2; j<Nj+2; ++j) {
      for (int k=2; k<Nk+2; ++k) {


	E[IND(i,j,k,1)] += dt * dtE[IND(i,j,k,1)];
	E[IND(i,j,k,2)] += dt * dtE[IND(i,j,k,2)];
	E[IND(i,j,k,3)] += dt * dtE[IND(i,j,k,3)];

	B[IND(i,j,k,1)] += dt * dtB[IND(i,j,k,1)];
	B[IND(i,j,k,2)] += dt * dtB[IND(i,j,k,2)];
	B[IND(i,j,k,3)] += dt * dtB[IND(i,j,k,3)];
      }
    }
  }



  sim->status.time_simulation += dt;
  sim->status.iteration += 1;
}



int main(int argc, char **argv)
{
  cow_init(0, NULL, 0);

  struct ffe_sim sim;

  ffe_sim_init(&sim);
  ffe_sim_initial_data(&sim);

  cow_dfield_write(sim.magnetic[0], "chkpt.0000.h5");
  cow_dfield_write(sim.electric[0], "chkpt.0000.h5");

  while (sim.status.time_simulation < 0.00001) {
    ffe_sim_advance(&sim);
    printf("n=%06d t=%3.2e\n", sim.status.iteration, sim.status.time_simulation);
  }

  cow_dfield_write(sim.magnetic[1], "chkpt.0001.h5");
  cow_dfield_write(sim.electric[1], "chkpt.0001.h5");


  ffe_sim_free(&sim);

  cow_finalize();
  return 0;
}
