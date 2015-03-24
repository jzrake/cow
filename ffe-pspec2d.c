#include <stdio.h>
#include "cow.h"


static void process_ffe_hdf5_file(char *filename);


int main(int argc, char **argv)
{

  if (argc == 1) {
    printf("usage: ffe-pspsec2d input.h5\n");
    return 0;
  }


  //cow_init(argc, argv, COW_NOREOPEN_STDOUT);
  cow_init(argc, argv, 0);

  for (int n=1; n<argc; ++n) {
    process_ffe_hdf5_file(argv[n]);
  }

  cow_finalize();  

  return 0;
}






void process_ffe_hdf5_file(char *filename)
{
  cow_domain *domain = cow_domain_new();
  cow_dfield *magnetic = cow_dfield_new();
  cow_dfield *electric = cow_dfield_new();
  cow_dfield *vecpoten = cow_dfield_new();
  cow_dfield *helicity = cow_dfield_new();

  cow_domain_readsize(domain, filename, "cell_primitive/B1");
  cow_domain_commit(domain);


  cow_dfield_setdomain(magnetic, domain);
  cow_dfield_setname(magnetic, "cell_primitive");
  cow_dfield_addmember(magnetic, "B1");
  cow_dfield_addmember(magnetic, "B2");
  cow_dfield_addmember(magnetic, "B3");
  cow_dfield_commit(magnetic);


  cow_dfield_setdomain(electric, domain);
  cow_dfield_setname(electric, "cell_primitive");
  cow_dfield_addmember(electric, "E1");
  cow_dfield_addmember(electric, "E2");
  cow_dfield_addmember(electric, "E3");
  cow_dfield_commit(electric);


  cow_dfield_setdomain(vecpoten, domain);
  cow_dfield_setname(vecpoten, "vector_potential");
  cow_dfield_addmember(vecpoten, "A1");
  cow_dfield_addmember(vecpoten, "A2");
  cow_dfield_addmember(vecpoten, "A3");
  cow_dfield_commit(vecpoten);


  cow_dfield_setdomain(helicity, domain);
  cow_dfield_setname(helicity, "helicity");
  cow_dfield_addmember(helicity, "H");
  cow_dfield_commit(helicity);


  cow_dfield_read(magnetic, filename);
  cow_dfield_read(electric, filename);
  cow_fft_inversecurl(magnetic, vecpoten);

  double *A = (double*) cow_dfield_getdatabuffer(magnetic);
  double *B = (double*) cow_dfield_getdatabuffer(vecpoten);
  double *H = (double*) cow_dfield_getdatabuffer(helicity);

  for (int n=0; n<cow_domain_getnumlocalzonesincguard(domain, COW_ALL_DIMS); ++n) {
    H[n] = A[3*n+0] * B[3*n+0] + A[3*n+1] * B[3*n+1] + A[3*n+2] * B[3*n+2];
  }



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


  cow_histogram *Hm = cow_histogram_new();
  cow_histogram_setlower(Hm, 0, 1);
  cow_histogram_setupper(Hm, 0, 8192);
  cow_histogram_setnbins(Hm, 0, 4096);
  cow_histogram_setspacing(Hm, COW_HIST_SPACING_LINEAR);
  cow_histogram_setfullname(Hm, "helicity");
  cow_histogram_setnickname(Hm, "helicity");



  cow_fft_pspecvecfield(magnetic, Pb);
  cow_fft_pspecvecfield(electric, Pe);
  cow_fft_pspecscafield(helicity, Hm);

  cow_histogram_dumphdf5(Pb, filename, "spectra");
  cow_histogram_dumphdf5(Pe, filename, "spectra");
  cow_histogram_dumphdf5(Hm, filename, "spectra");

  /* cow_histogram_dumpascii(Pb, "magnetic.dat"); */
  /* cow_histogram_dumpascii(Pe, "electric.dat"); */
  /* cow_histogram_dumpascii(Hm, "helicity.dat"); */

  cow_histogram_del(Pb);
  cow_histogram_del(Pe);
  cow_histogram_del(Hm);
  cow_dfield_del(magnetic);
  cow_dfield_del(electric);
  cow_dfield_del(vecpoten);
  cow_dfield_del(helicity);
  cow_domain_del(domain);
}
