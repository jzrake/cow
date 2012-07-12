
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "cow.h"
#if (COW_MPI)
#include <mpi.h>
#endif
#define KILOBYTES (1<<10)
#define MEGABYTES (1<<20)
#define GETENVINT(a,dflt) (getenv(a) ? atoi(getenv(a)) : dflt)
#define GETENVDBL(a,dflt) (getenv(a) ? atof(getenv(a)) : dflt)



static void div5(double *result, double **args, int **s, void *u)
{
#define diff5(f,s) ((-f[2*s] + 8*f[s] - 8*f[-s] + f[-2*s]) / 12.0)
  double *f0 = &args[0][0];
  double *f1 = &args[0][1];
  double *f2 = &args[0][2];
  *result = diff5(f0, s[0][0]) + diff5(f1, s[0][1]) + diff5(f2, s[0][2]);
#undef diff5
}
static void rot5(double *result, double **args, int **s, void *u)
{
  // http://en.wikipedia.org/wiki/Five-point_stencil
#define diff5(f,s) ((-f[2*s] + 8*f[s] - 8*f[-s] + f[-2*s]) / 12.0)
  double *f0 = &args[0][0];
  double *f1 = &args[0][1];
  double *f2 = &args[0][2];
  result[0] = diff5(f2, s[0][1]) - diff5(f1, s[0][2]);
  result[1] = diff5(f0, s[0][2]) - diff5(f2, s[0][0]);
  result[2] = diff5(f1, s[0][0]) - diff5(f0, s[0][1]);
#undef diff5
}
static void take_elm0(double *result, double **args, int **s, void *u)
{
  *result = args[0][0];
}
static void take_mag3(double *result, double **args, int **s, void *u)
{
  double *m = args[0];
  *result = sqrt(m[0]*m[0] + m[1]*m[1] + m[2]*m[2]);
}
static void take_lorentzfactor(double *result, double **args, int **s, void *u)
{
  double *m = args[0];
  *result = 1.0 / sqrt(1.0 - (m[0]*m[0] + m[1]*m[1] + m[2]*m[2]));
}

cow_dfield *cow_vectorfield(cow_domain *domain, const char *name)
{
  cow_dfield *f = cow_dfield_new(domain, name);
  cow_dfield_addmember(f, "0");
  cow_dfield_addmember(f, "1");
  cow_dfield_addmember(f, "2");
  cow_dfield_commit(f);
  return f;
}
cow_dfield *cow_scalarfield(cow_domain *domain, const char *name)
{
  cow_dfield *f = cow_dfield_new(domain, name);
  cow_dfield_addmember(f, "0");
  cow_dfield_commit(f);
  return f;
}

void make_hist(cow_dfield *f, cow_transform op, const char *fout, const char *m)
{
  char nickname[1024];
  snprintf(nickname, 1024, "%s-hist", m ? m : cow_dfield_getname(f));

  double reduc[3]; // min, max, sum
  cow_dfield_reduce(f, op, reduc);

  cow_histogram *hist = cow_histogram_new();
  cow_histogram_setlower(hist, 0, reduc[0]);
  cow_histogram_setupper(hist, 0, reduc[1]);
  cow_histogram_setnbins(hist, 0, 500);
  cow_histogram_setbinmode(hist, COW_HIST_BINMODE_COUNTS);
  cow_histogram_setdomaincomm(hist, cow_dfield_getdomain(f));
  cow_histogram_commit(hist);
  cow_histogram_setnickname(hist, nickname);
  cow_histogram_populate(hist, f, op);
  cow_histogram_dumphdf5(hist, fout, "");
  cow_histogram_del(hist);
}

int main(int argc, char **argv)
{
  int modes = 0;
  int collective = GETENVINT("COW_HDF5_COLLECTIVE", 0);
  modes |= GETENVINT("COW_NOREOPEN_STDOUT", 0) ? COW_NOREOPEN_STDOUT : 0;
  modes |= GETENVINT("COW_DISABLE_MPI", 0) ? COW_DISABLE_MPI : 0;

  cow_init(argc, argv, modes);
  if (argc == 3) {
    printf("running on input file %s\n", argv[1]);
  }
  else {
    printf("usage: $> srhdhist infile.h5 outfile.h5\n");
    goto done;
  }
  printf("COW_HDF5_COLLECTIVE: %d\n", collective);

  char *finp = argv[1];
  char *fout = argv[2];
  cow_domain *domain = cow_domain_new();
  cow_domain_readsize(domain, finp, "prim/vx");
  cow_domain_setguard(domain, 2);
  cow_domain_commit(domain);
  cow_domain_setchunk(domain, 1);
  cow_domain_setcollective(domain, collective);
  cow_domain_setalign(domain, 4*KILOBYTES, 4*MEGABYTES);

  cow_dfield *vel = cow_dfield_new(domain, "prim");
  cow_dfield *rho = cow_dfield_new(domain, "prim");
  cow_dfield_addmember(vel, "vx");
  cow_dfield_addmember(vel, "vy");
  cow_dfield_addmember(vel, "vz");
  cow_dfield_addmember(rho, "rho");
  cow_dfield_commit(vel);
  cow_dfield_commit(rho);
  cow_dfield_read(vel, finp);
  cow_dfield_read(rho, finp);

  make_hist(vel, take_lorentzfactor, fout, "gamma");
  make_hist(rho, take_elm0, fout, "rho");

  cow_dfield *divV = cow_scalarfield(domain, "divV");
  cow_dfield *rotV = cow_vectorfield(domain, "rotV");
  cow_dfield_transform(divV, &vel, 1, div5, NULL);
  cow_dfield_transform(rotV, &vel, 1, rot5, NULL);

  make_hist(divV, take_elm0, fout, NULL);
  make_hist(rotV, take_mag3, fout, NULL);

  cow_dfield_del(vel);
  cow_dfield_del(rho);
  cow_dfield_del(divV);
  cow_dfield_del(rotV);
  cow_domain_del(domain);

 done:
  cow_finalize();
  return 0;
}
