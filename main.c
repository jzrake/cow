

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if (COW_MPI)
#include <mpi.h>
#endif


// -----------------------------------------------------------------------------
//
// These prototypes constitute the C.O.W. interface
//
// -----------------------------------------------------------------------------
struct cow_dfield; // forward declarations (for opaque data structure)
struct cow_domain;
typedef struct cow_dfield cow_dfield;
typedef struct cow_domain cow_domain;

cow_domain *cow_domain_new();
void cow_domain_commit(cow_domain *d);
void cow_domain_del(cow_domain *d);
cow_dfield *cow_domain_addfield(cow_domain *d, const char *name);
cow_dfield *cow_domain_getfield(cow_domain *d, const char *name);
void cow_domain_setsize(cow_domain *d, int dim, int size);
void cow_domain_setndim(cow_domain *d, int ndim);
void cow_domain_setguard(cow_domain *d, int guard);
void cow_domain_setprocsizes(cow_domain *d, int dim, int size);
cow_dfield *cow_domain_iteratefields(cow_domain *d);
cow_dfield *cow_domain_nextfield(cow_domain *d);

cow_dfield *cow_dfield_new();
void cow_dfield_commit(cow_dfield *f);
void cow_dfield_del(cow_dfield *f);
void cow_dfield_addmember(cow_dfield *f, const char *name);
void cow_dfield_setname(cow_dfield *f, const char *name);
const char *cow_dfield_getname(cow_dfield *f);
const char *cow_dfield_iteratemembers(cow_dfield *f);
const char *cow_dfield_nextmember(cow_dfield *f);
int cow_domain_getnumlocalzones(cow_domain *d);


// -----------------------------------------------------------------------------
//
// cow_domain interface functions
//
// -----------------------------------------------------------------------------
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
  int n_fields; // number of data fields (dynamically adjustable)
  int field_iter; // index into data fields array used for iterating over them
  int committed; // true after cow_domain_commit called, locks out most changes
  cow_dfield **fields; // array of pointers to data fields

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
  MPI_Datatype *send_type; // chunk of data to be sent to respective neighbor
  MPI_Datatype *recv_type; // " "                 received from " "
#endif
} ;

struct cow_domain *cow_domain_new()
{
  cow_domain *d = (cow_domain*) malloc(sizeof(cow_domain));
  cow_domain dom = {
    .glb_lower = { 0.0, 0.0, 0.0 },
    .glb_upper = { 1.0, 1.0, 1.0 },
    .loc_lower = { 0.0, 0.0, 0.0 },
    .loc_upper = { 1.0, 1.0, 1.0 },
    .L_nint = { 1, 1, 1 },
    .L_ntot = { 1, 1, 1 },
    .L_strt = { 1, 1, 1 },
    .G_ntot = { 1, 1, 1 },
    .G_strt = { 1, 1, 1 },
    .n_dims = 1,
    .n_ghst = 1,
    .n_fields = 0,
    .field_iter = 0,
    .committed = 0,
    .fields = NULL,
#if (COW_MPI)
    .proc_sizes = { 0, 0, 0 },
    .proc_index = { 0, 0, 0 },
#endif
  } ;
  *d = dom;
  return d;
}
void cow_domain_del(cow_domain *d)
{
#if (COW_MPI)
  MPI_Comm_free(&d->mpi_cart);
#endif
  for (int n=0; n<d->n_fields; ++n) cow_dfield_del(d->fields[n]);
  free(d->fields);
  free(d);
}
cow_dfield *cow_domain_addfield(cow_domain *d, const char *name)
{
  d->n_fields++;
  d->fields = (cow_dfield**) realloc(d->fields,
                                     d->n_fields*sizeof(cow_dfield*));
  cow_dfield *f = cow_dfield_new(d);
  cow_dfield_setname(f, name);
  if (d->committed) {
    // -------------------------------------------------------------------------
    // If the domain has already been committed, then also commit the new data
    // field. This way new fields may be added and removed dynamically to an
    // already committed domain.
    // -------------------------------------------------------------------------
    cow_dfield_commit(f);
  }
  return d->fields[d->n_fields-1] = f;
}
cow_dfield *cow_domain_getfield(cow_domain *d, const char *name)
{
  for (int n=0; n<d->n_fields; ++n) {
    if (strcmp(cow_dfield_getname(d->fields[n]), name) == 0) {
      return d->fields[n];
    }
  }
  return NULL;
}
cow_dfield *cow_domain_iteratefields(cow_domain *d)
{
  d->field_iter = 0;
  return cow_domain_nextfield(d);
}
cow_dfield *cow_domain_nextfield(cow_domain *d)
{
  return d->field_iter++ < d->n_fields ? d->fields[d->field_iter-1] : NULL;
}
void cow_domain_setsize(cow_domain *d, int dim, int size)
{
  if (dim > 3 || d->committed) return;
  d->G_ntot[dim] = size;
}
void cow_domain_setndim(cow_domain *d, int ndim)
{
  if (ndim > 3 || d->committed) return;
  d->n_dims = ndim;
}
void cow_domain_setguard(cow_domain *d, int guard)
{
  if (guard < 0 || d->committed) return;
  d->n_ghst = guard;
}
void cow_domain_setprocsizes(cow_domain *d, int dim, int size)
{
#if (COW_MPI)
  if (dim > 3 || d->committed) return;
  d->proc_sizes[dim] = size;
#endif
}
void cow_domain_commit(cow_domain *d)
{
  if (d->committed) return;

#if (COW_MPI)
  int w[3] = { 1, 1, 1 }; // `wrap`, periodic in all directions
  int r = 1; // `reorder` allow MPI to choose a cart_rank != comm_rank

  MPI_Comm_rank(MPI_COMM_WORLD, &d->comm_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &d->comm_size);
  MPI_Dims_create(d->comm_size, d->n_dims, d->proc_sizes);
  MPI_Cart_create(MPI_COMM_WORLD, d->n_dims, d->proc_sizes, w, r, &d->mpi_cart);
  MPI_Comm_rank(d->mpi_cart, &d->cart_rank);
  MPI_Cart_coords(d->mpi_cart, d->cart_rank, d->n_dims, d->proc_index);

  for (int i=0; i<d->n_dims; ++i) {
    // -------------------------------------------------------------------------
    // The number of subgrid zones for dimension i needs to be non-uniform if
    // proc_sizes[i] does not divide the G_ntot[i]. For each dimension, we add a
    // zone to the first R subgrids , where R is given below:
    // -------------------------------------------------------------------------
    const int R = d->G_ntot[i] % d->proc_sizes[i];
    const int normal_size = d->G_ntot[i] / d->proc_sizes[i];
    const int augmnt_size = normal_size + 1;
    const int thisdm_size = (d->proc_index[i]<R) ? augmnt_size : normal_size;
    const double dx = (d->glb_upper[i] - d->glb_lower[i]) / d->G_ntot[i];

    d->L_nint[i] = thisdm_size;
    d->G_strt[i] = 0;
    for (int j=0; j<d->proc_index[i]; ++j) {
      d->G_strt[i] += (j<R) ? augmnt_size : normal_size;
    }
    d->loc_lower[i] = d->glb_lower[i] + dx *  d->G_strt[i];
    d->loc_upper[i] = d->glb_upper[i] + dx * (d->G_strt[i] + thisdm_size);
    d->L_ntot[i] = d->L_nint[i] + 2 * d->n_ghst;
    d->L_strt[i] = d->n_ghst;
  }

  printf("[cow] subgrid layout is (%d %d %d)\n",
         d->proc_sizes[0], d->proc_sizes[1], d->proc_sizes[2]);
  //  printf("[cow] this domain has size (%d %d %d)\n",
  //         d->L_nint[0], d->L_nint[1], d->L_nint[2]);
#else
  for (int i=0; i<d->n_dims; ++i) {
    d->L_nint[i] = d->G_ntot[i];
    d->L_ntot[i] = d->G_ntot[i] + 2 * d->n_ghst;
    d->L_strt[i] = d->n_ghst;
    d->G_strt[i] = 0;

    d->loc_lower[i] = d->glb_lower[i];
    d->loc_upper[i] = d->glb_upper[i];
  }
#endif

  for (int n=0; n<d->n_fields; ++n) {
    cow_dfield_commit(d->fields[n]);
  }
  d->committed = 1;
}
int cow_domain_getnumlocalzones(cow_domain *d)
{
  return d->L_ntot[0] * d->L_ntot[1] * d->L_ntot[2];
}

// -----------------------------------------------------------------------------
//
// cow_dfield interface functions
//
// -----------------------------------------------------------------------------
struct cow_dfield
{
  char *name;
  char **members;
  int member_iter;
  int n_members;
  void *data;
  int committed;
  cow_domain *domain;
} ;
cow_dfield *cow_dfield_new(cow_domain *domain)
{
  cow_dfield *f = (cow_dfield*) malloc(sizeof(cow_dfield));
  f->name = NULL;
  f->members = NULL;
  f->n_members = 0;
  f->member_iter = 0;
  f->data = NULL;
  f->committed = 0;
  f->domain = domain;
  return f;
}
void cow_dfield_del(cow_dfield *f)
{
  for (int n=0; n<f->n_members; ++n) free(f->members[n]);
  free(f->members);
  free(f->name);
  free(f->data);
  free(f);
}
void cow_dfield_setname(cow_dfield *f, const char *name)
{
  if (f->committed) return;
  f->name = (char*) realloc(f->name, strlen(name)+1);
  strcpy(f->name, name);
}
const char *cow_dfield_getname(cow_dfield *f)
{
  return f->name;
}
void cow_dfield_add_member(cow_dfield *f, const char *name)
{
  if (f->committed) return;
  f->n_members++;
  f->members = (char**) realloc(f->members, f->n_members*sizeof(char*));
  f->members[f->n_members-1] = (char*) malloc(strlen(name)+1);
  strcpy(f->members[f->n_members-1], name);
}
const char *cow_dfield_iteratemembers(cow_dfield *f)
{
  f->member_iter = 0;
  return cow_dfield_nextmember(f);
}
const char *cow_dfield_nextmember(cow_dfield *f)
{
  return f->member_iter++ < f->n_members ? f->members[f->member_iter-1] : NULL;
}
void cow_dfield_commit(cow_dfield *f)
{
  if (f->committed) return;
  const int n_zones = cow_domain_getnumlocalzones(f->domain);
  f->data = malloc(n_zones * f->n_members * sizeof(double));
  f->committed = 1;
}

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

  cow_dfield_add_member(prim, "vx");
  cow_dfield_add_member(prim, "vy");
  cow_dfield_add_member(prim, "vz");

  cow_dfield_add_member(magf, "Bx");
  cow_dfield_add_member(magf, "By");
  cow_dfield_add_member(magf, "Bz");

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

  cow_domain_del(domain);
#if (COW_MPI)
  MPI_Finalize();
#endif
  return 0;
}
