
#include <stdio.h>
#include "cow.h"
#if (COW_MPI)
#include <mpi.h>
#endif
#define KILOBYTES (1<<10)
#define MEGABYTES (1<<20)


static void divergence(double *result, double **args, int **s, cow_domain *d)
{
#define M(i,j,k) ((i)*si + (j)*sj + (k)*sk)
  double *fx = &args[0][0];
  double *fy = &args[0][1];
  double *fz = &args[0][2];
  int si = s[0][0];
  int sj = s[0][1];
  int sk = s[0][2];
  int i = 0;
  int j = 0;
  int k = 0;
  *result = ((fx[ M(i+1,j  ,k  ) ] + fx[ M(i+1,j+1,k  ) ] +
	      fx[ M(i+1,j  ,k+1) ] + fx[ M(i+1,j+1,k+1) ]) -
	     (fx[ M(i  ,j  ,k  ) ] + fx[ M(i  ,j+1,k  ) ] +
	      fx[ M(i  ,j  ,k+1) ] + fx[ M(i  ,j+1,k+1) ])) / 4.0
    +       ((fy[ M(i  ,j+1,k  ) ] + fy[ M(i  ,j+1,k+1) ] +
	      fy[ M(i+1,j+1,k  ) ] + fy[ M(i+1,j+1,k+1) ]) -
	     (fy[ M(i  ,j  ,k  ) ] + fy[ M(i  ,j  ,k+1) ] +
	      fy[ M(i+1,j  ,k  ) ] + fy[ M(i+1,j  ,k+1) ])) / 4.0
    +       ((fz[ M(i  ,j  ,k+1) ] + fz[ M(i+1,j  ,k+1) ] +
	      fz[ M(i  ,j+1,k+1) ] + fz[ M(i+1,j+1,k+1) ]) -
	     (fz[ M(i  ,j  ,k  ) ] + fz[ M(i+1,j  ,k  ) ] +
	      fz[ M(i  ,j+1,k  ) ] + fz[ M(i+1,j+1,k  ) ])) / 4.0;
}


int main(int argc, char **argv)
{
#if (COW_MPI)
  {
    int rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank != 0) freopen("/dev/null", "w", stdout);
    printf("was compiled with MPI support\n");
  }
#endif
  cow_domain *domain = cow_domain_new();
  cow_dfield *vel = cow_dfield_new(domain, "prim");
  cow_dfield *mag = cow_dfield_new(domain, "prim");
  cow_dfield *div = cow_dfield_new(domain, "divergence");

  cow_domain_setndim(domain, 3);
  cow_domain_setguard(domain, 2);
  cow_domain_setsize(domain, 0, 16);
  cow_domain_setsize(domain, 1, 16);
  cow_domain_setsize(domain, 2, 16);
  cow_domain_commit(domain);

  cow_domain_setchunk(domain, 1);
  cow_domain_setcollective(domain, 0);
  cow_domain_setalign(domain, 4*KILOBYTES, 4*MEGABYTES);

  cow_dfield_addmember(vel, "vx");
  cow_dfield_addmember(vel, "vy");
  cow_dfield_addmember(vel, "vz");
  cow_dfield_addmember(mag, "Bx");
  cow_dfield_addmember(mag, "By");
  cow_dfield_addmember(mag, "Bz");
  cow_dfield_addmember(div, "divB");

  cow_dfield_commit(vel);
  cow_dfield_commit(mag);
  cow_dfield_commit(div);

  if (argc > 1) {
    printf("running on input file %s\n", argv[1]);
  }
  cow_dfield_read(vel, argv[1]);
  cow_dfield_read(mag, argv[1]);

  cow_dfield_transform(div, &mag, 1, divergence);

  cow_dfield_write(div, "output.h5");

  cow_dfield_del(vel);
  cow_dfield_del(mag);
  cow_domain_del(domain);
#if (COW_MPI)
  MPI_Finalize();
#endif
  return 0;
}
