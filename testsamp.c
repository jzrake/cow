
#include <stdio.h>
#include "cow.h"
#if (COW_MPI)
#include <mpi.h>
#endif

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
  cow_dfield *data = cow_dfield_new(domain, "data");

  cow_domain_setndim(domain, 3);
  cow_domain_setguard(domain, 2);
  cow_domain_setsize(domain, 0, 4);
  cow_domain_setsize(domain, 1, 4);
  cow_domain_setsize(domain, 2, 4);
  cow_domain_commit(domain);

  cow_dfield_addmember(data, "d1");
  cow_dfield_addmember(data, "d2");
  cow_dfield_addmember(data, "d3");
  cow_dfield_commit(data);

  double *A = (double*) cow_dfield_getbuffer(data);
  for (int i=0; i<cow_domain_getnumlocalzonesincguard(domain, COW_ALL_DIMS);
       ++i) {
    A[3*i + 0] = 0.1;
    A[3*i + 1] = 0.2;
    A[3*i + 2] = 0.3;
  }
  cow_dfield_syncguard(data);

  double r[3] = {0.5, 0.5, 0.5};
  double sample[3];
  cow_dfield_sample(data, r, sample);
  printf("%f %f %f\n", sample[0], sample[1], sample[2]);

  cow_dfield_del(data);
  cow_domain_del(domain);

#if (COW_MPI)
  MPI_Finalize();
#endif
  return 0;
}
