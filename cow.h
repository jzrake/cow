

#ifndef COW_HEADER_INCLUDED
#define COW_HEADER_INCLUDED
#include <stdlib.h>

#ifdef COW_PRIVATE_DEFS
#if (COW_MPI)
#include <mpi.h>
#endif // COW_MPI
#if (COW_HDF5)
#include <hdf5.h>
#endif // COW_HDF5
#endif // COW_PRIVATE_DEFS


#define COW_ALL_DIMS -10
#define COW_HIST_SPACING_LINEAR -42
#define COW_HIST_SPACING_LOG -43
#define COW_HIST_BINMODE_DENSITY -44
#define COW_HIST_BINMODE_AVERAGE -45
#define COW_HIST_BINMODE_COUNTS -46


// -----------------------------------------------------------------------------
//
// These prototypes constitute the C.O.W. interface
//
// -----------------------------------------------------------------------------
struct cow_domain; // forward declarations (for opaque data structure)
struct cow_dfield;
typedef struct cow_domain cow_domain;
typedef struct cow_dfield cow_dfield;
typedef struct cow_histogram cow_histogram;
typedef void (*cow_transform)(double *result, double **args, int **strides,
			      void *udata);

void test_trans(double *result, double **args, int **strides, void *udata);


cow_domain *cow_domain_new();
void cow_domain_commit(cow_domain *d);
void cow_domain_del(cow_domain *d);
void cow_domain_setsize(cow_domain *d, int dim, int size);
void cow_domain_setndim(cow_domain *d, int ndim);
void cow_domain_setguard(cow_domain *d, int guard);
void cow_domain_setprocsizes(cow_domain *d, int dim, int size);
void cow_domain_setcollective(cow_domain *d, int mode);
void cow_domain_setchunk(cow_domain *d, int mode);
void cow_domain_setalign(cow_domain *d, int alignthreshold, int diskblocksize);
void cow_domain_readsize(cow_domain *d, const char *fname, const char *dname);
int cow_domain_getsize(cow_domain *d, int dim);
int cow_domain_getguard(cow_domain *d);
int cow_domain_getnumlocalzonesincguard(cow_domain *d, int dim);
int cow_domain_getnumlocalzonesinterior(cow_domain *d, int dim);
int cow_domain_getnumglobalzones(cow_domain *d, int dim);
int cow_domain_getglobalstartindex(cow_domain *d, int dim);

cow_dfield *cow_dfield_new(cow_domain *domain, const char *name);
void cow_dfield_commit(cow_dfield *f);
void cow_dfield_del(cow_dfield *f);
void cow_dfield_addmember(cow_dfield *f, const char *name);
void cow_dfield_setname(cow_dfield *f, const char *name);
void cow_dfield_extract(cow_dfield *f, const int *I0, const int *I1, void *out);
void cow_dfield_replace(cow_dfield *f, const int *I0, const int *I1, void *out);
void cow_dfield_loop(cow_dfield *f, cow_transform op, void *udata);
void cow_dfield_transform(cow_dfield *result, cow_dfield **args, int nargs,
			  cow_transform op, void *udata);

const char *cow_dfield_iteratemembers(cow_dfield *f);
const char *cow_dfield_nextmember(cow_dfield *f);
const char *cow_dfield_getname(cow_dfield *f);
int cow_dfield_getstride(cow_dfield *d, int dim);
size_t cow_dfield_getdatabytes(cow_dfield *f);
void *cow_dfield_getdata(cow_dfield *f);
void cow_dfield_syncguard(cow_dfield *f);
void cow_dfield_write(cow_dfield *f, const char *fname);
void cow_dfield_read(cow_dfield *f, const char *fname);

cow_histogram *cow_histogram_new();
void cow_histogram_commit(cow_histogram *h);
void cow_histogram_del(cow_histogram *h);
void cow_histogram_setbinmode(cow_histogram *h, int binmode);
void cow_histogram_setspacing(cow_histogram *h, int spacing);
void cow_histogram_setnbins(cow_histogram *h, int dim, int nbinsx);
void cow_histogram_setlower(cow_histogram *h, int dim, double v0);
void cow_histogram_setupper(cow_histogram *h, int dim, double v1);
void cow_histogram_setfullname(cow_histogram *h, const char *fullname);
void cow_histogram_setnickname(cow_histogram *h, const char *nickname);
void cow_histogram_setdomaincomm(cow_histogram *h, cow_domain *d);
void cow_histogram_addsample1(cow_histogram *h, double x, double w);
void cow_histogram_addsample2(cow_histogram *h, double x, double y, double w);
void cow_histogram_dumpascii(cow_histogram *h, const char *fn);
void cow_histogram_dumphdf5(cow_histogram *h, const char *fn, const char *dn);
void cow_histogram_synchronize(cow_histogram *h);
void cow_histogram_populate(cow_histogram *h, cow_dfield *f, cow_transform op);
double cow_histogram_getbinval(cow_histogram *h, int i, int j);



#ifdef COW_PRIVATE_DEFS

void _io_domain_commit(cow_domain *d);
void _io_domain_del(cow_domain *d);

struct cow_domain
{
  double glb_lower[3]; // lower coordinates of global physical domain
  double glb_upper[3]; // upper " "
  double loc_lower[3]; // lower coordinates of local physical domain
  double loc_upper[3]; // upper " "
  int L_nint[3]; // interior zones on the local subgrid
  int L_ntot[3]; // total " ", including guard zones
  int L_strt[3]; // starting index of interior zones on local subgrid
  int G_ntot[3]; // global domain size
  int G_strt[3]; // starting index into global domain
  int n_dims; // number of dimensions: 1, 2, 3
  int n_ghst; // number of guard zones: >= 0
  //  int n_fields; // number of data fields (dynamically adjustable)
  //  int field_iter; // index into data fields array used for iterating over them
  int balanced; // true when all subgrids have the same size
  int committed; // true after cow_domain_commit called, locks out size changes
  //  cow_dfield **fields; // array of pointers to data fields
#if (COW_MPI)
  int comm_rank; // rank with respect to MPI_COMM_WORLD communicator
  int comm_size; // size " "
  int cart_rank; // rank with respect to the cartesian communicator
  int cart_size; // size " "
  int proc_sizes[3]; // number of subgrids along each dimension
  int proc_index[3]; // coordinates of local subgrid in cartesian communicator
  int num_neighbors; // 3, 9, or 27 depending on the domain dimensionality
  int *neighbors; // cartesian ranks of the neighboring processors
  int *send_tags; // tag used to on send calls with respective neighbor
  int *recv_tags; // " "            recv " "
  MPI_Comm mpi_cart; // the cartesian communicator
#endif
#if (COW_HDF5)
  hsize_t L_nint_h5[3];
  hsize_t L_ntot_h5[3];
  hsize_t L_strt_h5[3];
  hsize_t G_ntot_h5[3];
  hsize_t G_strt_h5[3];
  hid_t fapl; // file access property list
  hid_t dcpl; // data set creation property list
  hid_t dxpl; // data set transfer property list
#endif
} ;

struct cow_dfield
{
  char *name;
  char **members;
  int member_iter;
  int n_members;
  void *data;
  int stride[3];
  int committed;
  cow_domain *domain;
#if (COW_MPI)
  MPI_Datatype *send_type; // chunk of data to be sent to respective neighbor
  MPI_Datatype *recv_type; // " "                 received from " "
#endif
} ;

struct cow_histogram
{
  int nbinsx;
  int nbinsy;
  double x0;
  double x1;
  double y0;
  double y1;
  double *bedgesx;
  double *bedgesy;
  double *weight;
  long *counts;
  char *nickname;
  char *fullname;
  int binmode;
  int spacing;
  int n_dims;
  int committed;
  cow_transform transform;
#if (COW_MPI)
  MPI_Comm comm;
#endif
} ;

#endif // COW_PRIVATE_DEFS
#endif // COW_HEADER_INCLUDED
