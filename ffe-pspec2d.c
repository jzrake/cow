#include <stdio.h>
#include "cow.h"



struct config_t
{
  char *dset_name_for_size;
  char *dset_name_primitive;
  char *dset_name_B1;
  char *dset_name_B2;
  char *dset_name_B3;
  char *dset_name_E1;
  char *dset_name_E2;
  char *dset_name_E3;
  int write_derived_fields;
} ;

static void process_file(char *filename, struct config_t cfg);


static struct config_t configure_for_ffe()
{
  struct config_t cfg;
  cfg.dset_name_for_size = NULL;
  cfg.dset_name_primitive = "cell_primitive";
  cfg.dset_name_B1 = "B1";
  cfg.dset_name_B2 = "B2";
  cfg.dset_name_B3 = "B3";
  cfg.dset_name_E1 = "E1";
  cfg.dset_name_E2 = "E2";
  cfg.dset_name_E3 = "E3";
  cfg.write_derived_fields = 0;
  return cfg;
}


static struct config_t configure_for_mara()
{
  struct config_t cfg;
  cfg.dset_name_for_size = "prim/Bx";
  cfg.dset_name_primitive = "prim";
  cfg.dset_name_B1 = "Bx";
  cfg.dset_name_B2 = "By";
  cfg.dset_name_B3 = "Bz";
  cfg.dset_name_E1 = "vx";
  cfg.dset_name_E2 = "vy";
  cfg.dset_name_E3 = "vz";
  cfg.write_derived_fields = 0;
  return cfg;
}


static struct config_t configure_for_m2()
{
  struct config_t cfg;
  cfg.dset_name_for_size = "cell_primitive/B1";
  cfg.dset_name_primitive = "cell_primitive";
  cfg.dset_name_B1 = "B1";
  cfg.dset_name_B2 = "B2";
  cfg.dset_name_B3 = "B3";
  cfg.dset_name_E1 = "v1";
  cfg.dset_name_E2 = "v2";
  cfg.dset_name_E3 = "v3";
  cfg.write_derived_fields = 0;
  return cfg;
}


enum DatasetType {
  DSET_TYPE_FFE,
  DSET_TYPE_MARA,
  DSET_TYPE_M2,
} ;


int main(int argc, char **argv)
{
  if (argc == 1) {
    printf("usage: ffe-pspsec2d input.h5\n");
    return 0;
  }


  enum DatasetType dset_type = DSET_TYPE_M2;
  struct config_t cfg;

  switch (dset_type) {
  case DSET_TYPE_FFE : cfg = configure_for_ffe();  break;
  case DSET_TYPE_MARA: cfg = configure_for_mara(); break;
  case DSET_TYPE_M2  : cfg = configure_for_m2();   break;
  default: cfg = configure_for_m2(); break;
  }

  cfg.write_derived_fields = 1;




  cow_init(0, NULL, COW_DISABLE_MPI);

  for (int n=1; n<argc; ++n) {
    process_file(argv[n], cfg);
  }

  cow_finalize();  

  return 0;
}






void process_file(char *filename, struct config_t cfg)
{
  cow_domain *domain = cow_domain_new();
  cow_dfield *magnetic = cow_dfield_new();
  cow_dfield *electric = cow_dfield_new();
  cow_dfield *vecpoten = cow_dfield_new();
  cow_dfield *jcurrent = cow_dfield_new();
  cow_dfield *helicity = cow_dfield_new();

  cow_domain_readsize(domain, filename, cfg.dset_name_for_size);
  cow_domain_commit(domain);


  cow_dfield_setdomain(magnetic, domain);
  cow_dfield_setname(magnetic, cfg.dset_name_primitive);
  cow_dfield_addmember(magnetic, cfg.dset_name_B1);
  cow_dfield_addmember(magnetic, cfg.dset_name_B2);
  cow_dfield_addmember(magnetic, cfg.dset_name_B3);
  cow_dfield_commit(magnetic);

  cow_dfield_setdomain(electric, domain);
  cow_dfield_setname(electric, cfg.dset_name_primitive);
  cow_dfield_addmember(electric, cfg.dset_name_E1);
  cow_dfield_addmember(electric, cfg.dset_name_E2);
  cow_dfield_addmember(electric, cfg.dset_name_E3);
  cow_dfield_commit(electric);

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

  cow_dfield_setdomain(helicity, domain);
  cow_dfield_setname(helicity, "helicity");
  cow_dfield_addmember(helicity, "H");
  cow_dfield_commit(helicity);



  cow_dfield_read(magnetic, filename);
  cow_dfield_read(electric, filename);
  cow_fft_inversecurl(magnetic, vecpoten);
  cow_fft_curl(magnetic, jcurrent);


  double *A = (double*) cow_dfield_getdatabuffer(vecpoten);
  double *B = (double*) cow_dfield_getdatabuffer(magnetic);
  double *H = (double*) cow_dfield_getdatabuffer(helicity);
  //double *E = (double*) cow_dfield_getdatabuffer(electric);
  //double *J = (double*) cow_dfield_getdatabuffer(jcurrent);


  double htot = 0.0;
  for (int n=0; n<cow_domain_getnumlocalzonesincguard(domain, COW_ALL_DIMS); ++n) {
    H[n] = A[3*n+0] * B[3*n+0] + A[3*n+1] * B[3*n+1] + A[3*n+2] * B[3*n+2];
    htot += H[n];
  }

  htot = cow_domain_dblsum(domain, htot) / cow_domain_getnumglobalzones(domain, COW_ALL_DIMS);
  printf("total magnetic helicity: %8.6e\n", htot);


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



  cow_fft_pspecvecfield(magnetic, Pb);
  cow_fft_pspecvecfield(electric, Pe);
  cow_fft_helicityspec(magnetic, Hr, Hi);

  cow_histogram_dumphdf5(Pb, filename, "spectra");
  cow_histogram_dumphdf5(Pe, filename, "spectra");
  cow_histogram_dumphdf5(Hi, filename, "spectra");
  cow_histogram_dumphdf5(Hr, filename, "spectra");

  if (cfg.write_derived_fields) {
    cow_dfield_write(vecpoten, filename);
    cow_dfield_write(jcurrent, filename);
  }

  cow_histogram_del(Pb);
  cow_histogram_del(Pe);
  cow_histogram_del(Hr);
  cow_histogram_del(Hi);
  cow_dfield_del(magnetic);
  cow_dfield_del(electric);
  cow_dfield_del(vecpoten);
  cow_dfield_del(jcurrent);
  cow_dfield_del(helicity);
  cow_domain_del(domain);
}
