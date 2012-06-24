
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
    printf("was compiled with COW_MPI\n");
  }
#endif

  cow_domain *domain = cow_domain_new();
  cow_dfield *prim = cow_domain_addfield(domain, "primitive");
  cow_dfield *magf = cow_domain_addfield(domain, "magnetic");

  cow_dfield_addmember(prim, "vx");
  cow_dfield_addmember(prim, "vy");
  cow_dfield_addmember(prim, "vz");

  cow_dfield_addmember(magf, "Bx");
  cow_dfield_addmember(magf, "By");
  cow_dfield_addmember(magf, "Bz");

  cow_domain_setndim(domain, 3);
  cow_domain_setguard(domain, 2);
  cow_domain_setsize(domain, 0, 12);
  cow_domain_setsize(domain, 1, 13);
  cow_domain_setsize(domain, 2, 14);
  cow_domain_commit(domain);

  for (cow_dfield *d = cow_domain_iteratefields(domain);
       d != NULL; d = cow_domain_nextfield(domain)) {
    printf("%s\n", cow_dfield_getname(d));
    for (const char *m = cow_dfield_iteratemembers(d);
         m != NULL; m = cow_dfield_nextmember(d)) {
      printf("\t%s\n", m);
    }
  }

  cow_dfield_syncguard(prim);
  cow_domain_del(domain);
#if (COW_MPI)
  MPI_Finalize();
#endif
  return 0;
}
