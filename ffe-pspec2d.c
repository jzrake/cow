#include <stdio.h>
#include <math.h>
#include <string.h>
#include "cow.h"
#include "read-ff-sdf.h"

#ifndef M_PI
#define M_PI 3.1415926535897926
#endif



enum DatasetType {
  DSET_TYPE_FFE=1,
  DSET_TYPE_MARA=2,
  DSET_TYPE_M2=3,
  DSET_TYPE_MFFE=4,
} ;

struct config_t
{
  enum DatasetType dset_type;
  char *filename;
  char *output_filename;
  char *dset_name_for_size;
  char *dset_name_electric;
  char *dset_name_magnetic;
  char *dset_name_B1;
  char *dset_name_B2;
  char *dset_name_B3;
  char *dset_name_E1;
  char *dset_name_E2;
  char *dset_name_E3;
  int write_derived_fields;
  int resolution; /* only used for sdf files; array size not inferred */
  int three_d;
  int num_blocks;
  int checkpoint_number;
  int checkpoint_number0;
  int checkpoint_number1;
  int max_pspec_bin;
  int num_pspec_bins;
} ;


static void process_file(struct config_t cfg);
static void read_data_from_file(struct config_t cfg, cow_dfield *magnetic, cow_dfield *electric, double *time);


static struct config_t configure_new()
{
  struct config_t cfg;
  cfg.dset_type = 0;
  cfg.filename = NULL;
  cfg.output_filename = NULL;
  cfg.dset_name_for_size = NULL;
  cfg.dset_name_electric = NULL;
  cfg.dset_name_magnetic = NULL;
  cfg.dset_name_B1 = NULL;
  cfg.dset_name_B2 = NULL;
  cfg.dset_name_B3 = NULL;
  cfg.dset_name_E1 = NULL;
  cfg.dset_name_E2 = NULL;
  cfg.dset_name_E3 = NULL;
  cfg.write_derived_fields = 0;
  cfg.resolution = 0;
  cfg.three_d = 0;
  cfg.num_blocks = 0;
  cfg.checkpoint_number = 0;
  cfg.checkpoint_number0 = 0;
  cfg.checkpoint_number1 = 128;
  cfg.max_pspec_bin = 8192;
  cfg.num_pspec_bins = 4096;
  return cfg;
}

static void configure_for_ffe(struct config_t *cfg)
{
  cfg->dset_type = DSET_TYPE_FFE;
  cfg->dset_name_electric = "cell_primitive";
  cfg->dset_name_magnetic = "cell_primitive";
  cfg->dset_name_B1 = "B1";
  cfg->dset_name_B2 = "B2";
  cfg->dset_name_B3 = "B3";
  cfg->dset_name_E1 = "E1";
  cfg->dset_name_E2 = "E2";
  cfg->dset_name_E3 = "E3";
}

static void configure_for_mara(struct config_t *cfg)
{
  cfg->dset_type = DSET_TYPE_MARA;
  cfg->dset_name_for_size = "prim/Bx";
  cfg->dset_name_electric = "prim";
  cfg->dset_name_magnetic = "prim";
  cfg->dset_name_B1 = "Bx";
  cfg->dset_name_B2 = "By";
  cfg->dset_name_B3 = "Bz";
  cfg->dset_name_E1 = "vx";
  cfg->dset_name_E2 = "vy";
  cfg->dset_name_E3 = "vz";
}

static void configure_for_m2(struct config_t *cfg)
{
  cfg->dset_type = DSET_TYPE_M2;
  cfg->dset_name_for_size = "cell_primitive/B1";
  cfg->dset_name_electric = "cell_primitive";
  cfg->dset_name_magnetic = "cell_primitive";
  cfg->dset_name_B1 = "B1";
  cfg->dset_name_B2 = "B2";
  cfg->dset_name_B3 = "B3";
  cfg->dset_name_E1 = "v1";
  cfg->dset_name_E2 = "v2";
  cfg->dset_name_E3 = "v3";
}

static void configure_for_mffe(struct config_t *cfg)
{
  cfg->dset_type = DSET_TYPE_MFFE;
  cfg->dset_name_for_size = "magnetic/B1";
  cfg->dset_name_electric = "electric";
  cfg->dset_name_magnetic = "magnetic";
  cfg->dset_name_B1 = "B1";
  cfg->dset_name_B2 = "B2";
  cfg->dset_name_B3 = "B3";
  cfg->dset_name_E1 = "E1";
  cfg->dset_name_E2 = "E2";
  cfg->dset_name_E3 = "E3";
}



int main(int argc, char **argv)
{
  if (argc == 1) {
    printf("usage: ffe-pspsec2d input.h5 [options]\n");
    return 0;
  }



  int modes = 0;

  for (int n=2; n<argc; ++n) {

    if (argv[n][0] == '-') {
      if (argv[n][1] == 's') {
	printf("[cfg] running in serial mode\n");
	modes |= COW_DISABLE_MPI;
      }
      else {
	printf("[cfg] error: unknown arugment '%s'\n", argv[n]);
	return 1;
      }
    }
  }


  cow_init(0, NULL, modes);


  enum DatasetType dset_type = DSET_TYPE_MFFE;
  struct config_t cfg = configure_new();
  char output_filename[1024] = "";

  for (int n=2; n<argc; ++n) {

    if (argv[n][0] == '-') {
      continue; /* this was an option, not an argument*/
    }
    else if (strchr(argv[n], '=') == NULL) {
      continue; /* this was a filename, not an argument */
    }

    else if (!strncmp(argv[n], "output=", 7)) {
      sscanf(argv[n], "output=%1024s", output_filename);
      cfg.output_filename = output_filename;
    }
    else if (!strcmp(argv[n], "format=ffe")) {
      dset_type = DSET_TYPE_FFE;
    }
    else if (!strcmp(argv[n], "format=m2")) {
      dset_type = DSET_TYPE_M2;
    }
    else if (!strcmp(argv[n], "format=mara")) {
      dset_type = DSET_TYPE_MARA;
    }
    else if (!strcmp(argv[n], "format=mffe")) {
      dset_type = DSET_TYPE_MFFE;
    }
    else if (!strncmp(argv[n], "resolution=", 11)) {
      sscanf(argv[n], "resolution=%d", &cfg.resolution);
    }
    else if (!strncmp(argv[n], "num_blocks=", 11)) {
      sscanf(argv[n], "num_blocks=%d", &cfg.num_blocks);
    }
    else if (!strncmp(argv[n], "checkpoint_numbers=", 19)) {
      sscanf(argv[n], "checkpoint_numbers=%d-%d",
	     &cfg.checkpoint_number0, &cfg.checkpoint_number1);
    }
    else if (!strncmp(argv[n], "write_derived_fields=", 21)) {
      sscanf(argv[n], "write_derived_fields=%d", &cfg.write_derived_fields);
    }
    else {
      printf("[cfg] error: unknown option '%s'\n", argv[n]);
      return 1;
    }


  }


  printf("[cfg] using output=%s\n", cfg.output_filename);
  printf("[cfg] using resolution=%d\n", cfg.resolution);
  printf("[cfg] using three_d=%d\n", cfg.three_d);
  printf("[cfg] using num_blocks=%d\n", cfg.num_blocks);
  printf("[cfg] using write_derived_fields=%d\n", cfg.write_derived_fields);
  printf("[cfg] using num_pspec_bins=%d\n", cfg.num_pspec_bins);
  printf("[cfg] using max_pspec_bin=%d\n", cfg.max_pspec_bin);
  printf("[cfg] using checkpoint_numbers=%d-%d\n",
	 cfg.checkpoint_number0, cfg.checkpoint_number1);


  /* ================================================================
   * ================= setup for input data format ==================
   * ================================================================ */
  switch (dset_type) {
  case DSET_TYPE_FFE:
    configure_for_ffe(&cfg);
    printf("[cfg] using format=ffe\n");
    break;
  case DSET_TYPE_MARA:
    configure_for_mara(&cfg);
    printf("[cfg] using format=mara\n");
    break;
  case DSET_TYPE_M2:
    configure_for_m2(&cfg);
    printf("[cfg] using format=m2\n");
    break;
  case DSET_TYPE_MFFE:
    configure_for_mffe(&cfg);
    printf("[cfg] using format=mffe\n");
    break;
  }


  /* ================================================================
   * ================= process input files ==========================
   * ================================================================ */

  if (cfg.dset_type != DSET_TYPE_FFE) {

    for (int n=1; n<argc; ++n) {

      if (strchr(argv[n], '=') != NULL) {
	continue; /* this was an argument, not a filename */
      }
      if (argv[n][0] == '-') {
	continue; /* this was an option, not a filename */
      }

      cfg.filename = argv[n];

      //if (cfg.output_filename == NULL) {
      cfg.output_filename = argv[n];
      //}

      process_file(cfg);

    }
  }

  else {

    for (int n=cfg.checkpoint_number0; n<=cfg.checkpoint_number1; ++n) {

      cfg.filename = argv[1];
      cfg.checkpoint_number = n;
      process_file(cfg);

    }
  }

  cow_finalize();

  return 0;
}






void process_file(struct config_t cfg)
{
  cow_domain *domain = cow_domain_new();
  cow_dfield *magnetic = cow_dfield_new();
  cow_dfield *electric = cow_dfield_new();
  cow_dfield *vecpoten = cow_dfield_new();
  cow_dfield *jcurrent = cow_dfield_new();
  cow_dfield *helicity = cow_dfield_new();


  /* Domain setup */
  /* ---------------------------------------------------- */
  if (cfg.dset_name_for_size) {
    cow_domain_readsize(domain, cfg.filename, cfg.dset_name_for_size);
  }
  else if (cfg.three_d) {
    cow_domain_setndim(domain, 3);
    cow_domain_setsize(domain, 0, cfg.resolution);
    cow_domain_setsize(domain, 1, cfg.resolution);
    cow_domain_setsize(domain, 2, cfg.resolution);
  }
  else {
    cow_domain_setndim(domain, 2);
    cow_domain_setsize(domain, 0, cfg.resolution);
    cow_domain_setsize(domain, 1, cfg.resolution);
  }
  cow_domain_commit(domain);


  /* Data fields setup */
  /* ---------------------------------------------------- */
  cow_dfield_setdomain(magnetic, domain);
  cow_dfield_setname(magnetic, cfg.dset_name_magnetic);
  cow_dfield_addmember(magnetic, cfg.dset_name_B1);
  cow_dfield_addmember(magnetic, cfg.dset_name_B2);
  cow_dfield_addmember(magnetic, cfg.dset_name_B3);
  cow_dfield_commit(magnetic);

  cow_dfield_setdomain(electric, domain);
  cow_dfield_setname(electric, cfg.dset_name_electric);
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


  /* Histograms setup */
  /* ---------------------------------------------------- */
  cow_histogram *Pb = cow_histogram_new();
  cow_histogram_setlower(Pb, 0, 1);
  cow_histogram_setupper(Pb, 0, cfg.max_pspec_bin);
  cow_histogram_setnbins(Pb, 0, cfg.num_pspec_bins);
  cow_histogram_setspacing(Pb, COW_HIST_SPACING_LINEAR);
  cow_histogram_setfullname(Pb, "magnetic");
  cow_histogram_setnickname(Pb, "magnetic");


  cow_histogram *Pe = cow_histogram_new();
  cow_histogram_setlower(Pe, 0, 1);
  cow_histogram_setupper(Pe, 0, cfg.max_pspec_bin);
  cow_histogram_setnbins(Pe, 0, cfg.num_pspec_bins);
  cow_histogram_setspacing(Pe, COW_HIST_SPACING_LINEAR);
  cow_histogram_setfullname(Pe, "electric");
  cow_histogram_setnickname(Pe, "electric");


  cow_histogram *Hr = cow_histogram_new();
  cow_histogram_setlower(Hr, 0, 1);
  cow_histogram_setupper(Hr, 0, cfg.max_pspec_bin);
  cow_histogram_setnbins(Hr, 0, cfg.num_pspec_bins);
  cow_histogram_setspacing(Hr, COW_HIST_SPACING_LINEAR);
  cow_histogram_setfullname(Hr, "helicity");
  cow_histogram_setnickname(Hr, "helicity");


  cow_histogram *Hi = cow_histogram_new();
  cow_histogram_setlower(Hi, 0, 1);
  cow_histogram_setupper(Hi, 0, cfg.max_pspec_bin);
  cow_histogram_setnbins(Hi, 0, cfg.num_pspec_bins);
  cow_histogram_setspacing(Hi, COW_HIST_SPACING_LINEAR);
  cow_histogram_setfullname(Hi, "helicity-imag");
  cow_histogram_setnickname(Hi, "helicity-imag");


  cow_histogram *alpha_hist = cow_histogram_new();
  cow_histogram_setlower(alpha_hist, 0, -512.0);
  cow_histogram_setupper(alpha_hist, 0, +512.0);
  cow_histogram_setnbins(alpha_hist, 0, 8192);
  cow_histogram_setspacing(alpha_hist, COW_HIST_SPACING_LINEAR);
  cow_histogram_setfullname(alpha_hist, "alpha-hist");
  cow_histogram_setnickname(alpha_hist, "alpha-hist");
  cow_histogram_setdomaincomm(alpha_hist, domain);
  cow_histogram_commit(alpha_hist);

  cow_histogram *mup_hist = cow_histogram_new(); /* 1 - cos(theta) */
  cow_histogram_setlower(mup_hist, 0, 1e-8);
  cow_histogram_setupper(mup_hist, 0, 1.0);
  cow_histogram_setnbins(mup_hist, 0, 512);
  cow_histogram_setspacing(mup_hist, COW_HIST_SPACING_LOG);
  cow_histogram_setbinmode(mup_hist, COW_HIST_BINMODE_DENSITY);
  cow_histogram_setfullname(mup_hist, "mup-hist");
  cow_histogram_setnickname(mup_hist, "mup-hist");
  cow_histogram_setdomaincomm(mup_hist, domain);
  cow_histogram_commit(mup_hist);


  cow_histogram *mum_hist = cow_histogram_new(); /* 1 + cos(theta) */
  cow_histogram_setlower(mum_hist, 0, 1e-8);
  cow_histogram_setupper(mum_hist, 0, 1.0);
  cow_histogram_setnbins(mum_hist, 0, 512);
  cow_histogram_setspacing(mum_hist, COW_HIST_SPACING_LOG);
  cow_histogram_setbinmode(mum_hist, COW_HIST_BINMODE_DENSITY);
  cow_histogram_setfullname(mum_hist, "mum-hist");
  cow_histogram_setnickname(mum_hist, "mum-hist");
  cow_histogram_setdomaincomm(mum_hist, domain);
  cow_histogram_commit(mum_hist);


  double time = 0.0;
  read_data_from_file(cfg, magnetic, electric, &time);

  cow_fft_inversecurl(magnetic, vecpoten);
  cow_fft_curl(magnetic, jcurrent);


  double *A = (double*) cow_dfield_getdatabuffer(vecpoten);
  double *B = (double*) cow_dfield_getdatabuffer(magnetic);
  double *H = (double*) cow_dfield_getdatabuffer(helicity);
  //double *E = (double*) cow_dfield_getdatabuffer(electric);
  double *J = (double*) cow_dfield_getdatabuffer(jcurrent);


  double htot = 0.0;
  for (int n=0; n<cow_domain_getnumlocalzonesincguard(domain, COW_ALL_DIMS); ++n) {
    H[n] = A[3*n+0] * B[3*n+0] + A[3*n+1] * B[3*n+1] + A[3*n+2] * B[3*n+2];
    htot += H[n];
  }

  for (int n=0; n<cow_domain_getnumlocalzonesincguard(domain, COW_ALL_DIMS); ++n) {
    double BB = B[3*n+0] * B[3*n+0] + B[3*n+1] * B[3*n+1] + B[3*n+2] * B[3*n+2];
    double JJ = J[3*n+0] * J[3*n+0] + J[3*n+1] * J[3*n+1] + J[3*n+2] * J[3*n+2];
    double JB = J[3*n+0] * B[3*n+0] + J[3*n+1] * B[3*n+1] + J[3*n+2] * B[3*n+2];
    cow_histogram_addsample1(alpha_hist, JB/BB/(2*M_PI), 1.0);
    cow_histogram_addsample1(mup_hist, 1 - JB / sqrt(BB * JJ), 1.0);
    cow_histogram_addsample1(mum_hist, 1 + JB / sqrt(BB * JJ), 1.0);
  }
  cow_histogram_seal(alpha_hist);
  cow_histogram_seal(mup_hist);
  cow_histogram_seal(mum_hist);


  htot = cow_domain_dblsum(domain, htot) / cow_domain_getnumglobalzones(domain, COW_ALL_DIMS);
  printf("[main] total magnetic helicity: %8.6e\n", htot);


  //cow_fft_pspecvecfield(magnetic, Pb);
  //cow_fft_pspecvecfield(electric, Pe);
  //cow_fft_helicityspec(magnetic, Hr, Hi);

  if (cfg.output_filename) {

    char gname[1024];

    if (cfg.dset_type == DSET_TYPE_FFE) {
      snprintf(gname, 1024, "spectra-%6.4e", time);
    }
    else {
      snprintf(gname, 1024, "spectra");
    }

    cow_histogram_dumphdf5(Pb, cfg.output_filename, gname);
    cow_histogram_dumphdf5(Pe, cfg.output_filename, gname);
    cow_histogram_dumphdf5(Hr, cfg.output_filename, gname);
    cow_histogram_dumphdf5(alpha_hist, cfg.output_filename, gname);
    cow_histogram_dumphdf5(mup_hist, cfg.output_filename, gname);
    cow_histogram_dumphdf5(mum_hist, cfg.output_filename, gname);
  }

  if (cfg.write_derived_fields && cfg.output_filename) {
    cow_dfield_write(vecpoten, cfg.output_filename);
    cow_dfield_write(jcurrent, cfg.output_filename);
  }

  cow_histogram_del(Pb);
  cow_histogram_del(Pe);
  cow_histogram_del(Hr);
  cow_histogram_del(Hi);
  cow_histogram_del(alpha_hist);
  cow_histogram_del(mup_hist);
  cow_histogram_del(mum_hist);
  cow_dfield_del(magnetic);
  cow_dfield_del(electric);
  cow_dfield_del(vecpoten);
  cow_dfield_del(jcurrent);
  cow_dfield_del(helicity);
  cow_domain_del(domain);
}



void read_data_from_file(struct config_t cfg,
			 cow_dfield *magnetic,
			 cow_dfield *electric, double *time)
{
  cow_domain *domain = cow_dfield_getdomain(magnetic);

  /* If it's Mara, M2, or MFFE */
  /* ---------------------------------------------------- */
  if (cfg.dset_type == DSET_TYPE_MARA ||
      cfg.dset_type == DSET_TYPE_M2 ||
      cfg.dset_type == DSET_TYPE_MFFE) {

    cow_dfield_read(magnetic, cfg.filename);
    cow_dfield_read(electric, cfg.filename);
    *time = 0.0;

  }

  /* If it's FFE */
  /* ---------------------------------------------------- */
  else if (cfg.dset_type == DSET_TYPE_FFE) {


    printf("[main] reading ffe dataset %s\n", cfg.filename);

    int mpi_rank = cow_domain_getcartrank(domain);
    int mpi_size = cow_domain_getcartsize(domain);
    int block0 = (mpi_rank + 0) * cfg.num_blocks / mpi_size;
    int block1 = (mpi_rank + 1) * cfg.num_blocks / mpi_size;


    if (cfg.num_blocks % mpi_size != 0 ||
	cfg.num_blocks < mpi_size) {
      printf("[main] error: MPI size must divide, and be smaller than or equal "
	     "to the number of blocks\n");
      return;
    }

    for (int n=block0; n<block1; ++n) {

      double *EBdata[6];
      int startIndices[3] = {0, 0, 0};
      int sizes[3] = {1, 1, 1};


      int error = getEBFromSDF(n, cfg.checkpoint_number, cfg.filename,
			       EBdata, startIndices, sizes, time);

      if (cow_domain_intsum(domain, error)) {
	printf("[main] error: something went wrong reading an SDF file\n");
	return;
      }

      int num_zones = sizes[0] * sizes[1];
      double *Edata = (double*) malloc(num_zones*3*sizeof(double));
      double *Bdata = (double*) malloc(num_zones*3*sizeof(double));


      for (int i=0; i<sizes[0]; ++i) {
	for (int j=0; j<sizes[1]; ++j) {

	  int ms = j * sizes[0] + i; /* source index (F) */
	  int md = i * sizes[1] + j; /* dest index (C) */

	  Edata[3*md + 0] = EBdata[0][ms];
	  Edata[3*md + 1] = EBdata[1][ms];
	  Edata[3*md + 2] = EBdata[2][ms];
	  Bdata[3*md + 0] = EBdata[3][ms];
	  Bdata[3*md + 1] = EBdata[4][ms];
	  Bdata[3*md + 2] = EBdata[5][ms];

	}
      }

      cow_dfield_remap(electric, startIndices, sizes, Edata);
      cow_dfield_remap(magnetic, startIndices, sizes, Bdata);

      printf("[main] checkpoint %d (t=%3.2e) block %d: i=%d j=%d nx=%d ny=%d\n",
	     cfg.checkpoint_number,
	     *time, n, startIndices[0], startIndices[1], sizes[0], sizes[1]);


      free(EBdata[0]);
      free(EBdata[1]);
      free(EBdata[2]);
      free(EBdata[3]);
      free(EBdata[4]);
      free(EBdata[5]);
      free(Edata);
      free(Bdata);
    }
  }
}
