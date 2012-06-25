
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

  cow_domain_setndim(domain, 1);
  cow_domain_setguard(domain, 3);
  cow_domain_setsize(domain, 0, 10);
  cow_domain_commit(domain);

  for (cow_dfield *d = cow_domain_iteratefields(domain);
       d != NULL; d = cow_domain_nextfield(domain)) {
    printf("%s\n", cow_dfield_getname(d));
    for (const char *m = cow_dfield_iteratemembers(d);
         m != NULL; m = cow_dfield_nextmember(d)) {
      printf("\t%s\n", m);
    }
  }

  int si = cow_dfield_getstride(prim, 0);
  int ng = cow_domain_getguard(domain);

  double *P = (double*) cow_dfield_getdata(prim);
  for (int i=ng; i<cow_domain_getsize(domain, 0)+ng; ++i) {
    P[si*i + 0] = 1.0;
    P[si*i + 1] = 2.0;
    P[si*i + 2] = 3.0;
    printf("(%02d) %f %f %f\n", i, P[si*i + 0], P[si*i + 1], P[si*i + 2]);
  }

  cow_dfield_syncguard(prim);
  printf("\n\n ------------------------------------------------- \n\n");
  for (int i=0; i<cow_domain_getsize(domain, 0)+2*ng; ++i) {
    printf("(%02d) %f %f %f\n", i, P[si*i + 0], P[si*i + 1], P[si*i + 2]);
  }

  int I0[] = { 0, 0, 0 };
  int I1[] = { 2, 0, 0 };
  double *subarray = (double*) malloc(2 * 3 * sizeof(double));
  cow_dfield_extract(prim, I0, I1, subarray);
  printf("%f %f\n", subarray[0], subarray[3]);
  printf("%f %f\n", subarray[1], subarray[4]);
  printf("%f %f\n", subarray[2], subarray[5]);
  subarray[0] = 10.0;
  cow_dfield_replace(prim, I0, I1, subarray);
  printf("%f %f\n", subarray[0], subarray[1]);
  free(subarray);

  cow_domain_del(domain);
#if (COW_MPI)
  MPI_Finalize();
#endif
  return 0;
}
