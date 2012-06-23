

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if (COW_MPI)
#include <mpi.h>
#endif

struct cow_dfield; // forward declarations (for opaque data structure)
struct cow_domain;

struct cow_dfield *cow_domain_addfield(struct cow_domain *d, const char *name);
struct cow_dfield *cow_domain_getfield(struct cow_domain *d, const char *name);
struct cow_dfield *cow_dfield_new();
struct cow_domain *cow_domain_new();
void cow_dfield_del(struct cow_dfield *f);
void cow_dfield_addmember(struct cow_dfield *f, const char *name);
void cow_dfield_setname(struct cow_dfield *f, const char *name);
const char *cow_dfield_getname(struct cow_dfield *f);

void cow_domain_setsize(struct cow_domain *d, int dim, int size);
void cow_domain_setndim(struct cow_domain *d, int ndim);
void cow_domain_setguard(struct cow_domain *d, int guard);
void cow_domain_commit(struct cow_domain *d);
struct cow_dfield *cow_domain_iteratefields(struct cow_domain *d);
struct cow_dfield *cow_domain_nextfield(struct cow_domain *d);

const char *cow_dfield_iteratemembers(struct cow_dfield *f);
const char *cow_dfield_nextmember(struct cow_dfield *f);

void cow_dfield_commit(struct cow_dfield *f);
int cow_domain_getnumlocalzones(struct cow_domain *d);

// -----------------------------------------------------------------------------
//
// struct cow_domain interface functions
//
// -----------------------------------------------------------------------------
struct cow_domain
{
  double coord_lower[3]; // lower coordinates of physical domain
  double coord_upper[3]; // upper " "
  int A_nint[3];
  int L_ntot[3];
  int L_strt[3];
  int G_ntot[3];
  int G_strt[3];
  int n_dims;
  int n_ghst;
  int n_fields;
  int field_iter;
  int committed;
  struct cow_dfield **fields; // array of pointers to data fields

#if (COW_MPI)
  int comm_rank; // rank with respect to MPI_COMM_WORLD communicator
  int comm_size; // size " "
  int cart_rank; // rank with respect to the cartesian communicator
  int cart_size; // size " "
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
  struct cow_domain *d = (struct cow_domain*) malloc(sizeof(struct cow_domain));
  struct cow_domain dom = {
    .coord_lower = { 0.0, 0.0, 0.0 },
    .coord_upper = { 1.0, 1.0, 1.0 },
    .A_nint = { 1, 1, 1 },
    .L_ntot = { 1, 1, 1 },
    .L_strt = { 1, 1, 1 },
    .G_ntot = { 1, 1, 1 },
    .G_strt = { 1, 1, 1 },
    .n_dims = 1,
    .n_ghst = 1,
    .n_fields = 0,
    .field_iter = 0,
    .committed = 0,
    .fields = NULL
  } ;
  *d = dom;
  return d;
}
void cow_domain_del(struct cow_domain *d)
{
  for (int n=0; n<d->n_fields; ++n) cow_dfield_del(d->fields[n]);
  free(d->fields);
  free(d);
}
struct cow_dfield *cow_domain_addfield(struct cow_domain *d, const char *name)
{
  d->n_fields++;
  d->fields = (struct cow_dfield**) realloc(d->fields,
					    d->n_fields*sizeof(struct cow_dfield*));
  struct cow_dfield *f = cow_dfield_new(d);
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
struct cow_dfield *cow_domain_getfield(struct cow_domain *d, const char *name)
{
  for (int n=0; n<d->n_fields; ++n) {
    if (strcmp(cow_dfield_getname(d->fields[n]), name) == 0) {
      return d->fields[n];
    }
  }
  return NULL;
}
struct cow_dfield *cow_domain_iteratefields(struct cow_domain *d)
{
  d->field_iter = 0;
  return cow_domain_nextfield(d);
}
struct cow_dfield *cow_domain_nextfield(struct cow_domain *d)
{
  return d->field_iter++ < d->n_fields ? d->fields[d->field_iter-1] : NULL;
}
void cow_domain_setsize(struct cow_domain *d, int dim, int size)
{
  if (dim > 3) return;
  d->G_ntot[dim] = size;
}
void cow_domain_setndim(struct cow_domain *d, int ndim)
{
  if (ndim > 3) return;
  d->n_dims = ndim;
}
void cow_domain_setguard(struct cow_domain *d, int guard)
{
  d->n_ghst = guard;
}
void cow_domain_commit(struct cow_domain *d)
{
  if (d->committed) return;

#if (COW_MPI)
  MPI_Comm_rank(MPI_COMM_WORLD, &d->comm_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &d->comm_size);
#endif

  d->L_ntot[0] = d->G_ntot[0] + 2 * d->n_ghst;
  d->L_ntot[1] = d->G_ntot[1] + 2 * d->n_ghst;
  d->L_ntot[2] = d->G_ntot[2] + 2 * d->n_ghst;

  for (int n=0; n<d->n_fields; ++n) {
    cow_dfield_commit(d->fields[n]);
  }
  d->committed = 1;
}
int cow_domain_getnumlocalzones(struct cow_domain *d)
{
  return d->L_ntot[0] * d->L_ntot[1] * d->L_ntot[2];
}

// -----------------------------------------------------------------------------
//
// struct cow_dfield interface functions
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
  struct cow_domain *domain;
} ;
struct cow_dfield *cow_dfield_new(struct cow_domain *domain)
{
  struct cow_dfield *f = (struct cow_dfield*) malloc(sizeof(struct cow_dfield));
  f->name = NULL;
  f->members = NULL;
  f->n_members = 0;
  f->member_iter = 0;
  f->data = NULL;
  f->committed = 0;
  f->domain = domain;
  return f;
}
void cow_dfield_del(struct cow_dfield *f)
{
  for (int n=0; n<f->n_members; ++n) free(f->members[n]);
  free(f->members);
  free(f->name);
  free(f->data);
  free(f);
}
void cow_dfield_setname(struct cow_dfield *f, const char *name)
{
  if (f->committed) return;
  f->name = (char*) realloc(f->name, strlen(name)+1);
  strcpy(f->name, name);
}
const char *cow_dfield_getname(struct cow_dfield *f)
{
  return f->name;
}
void cow_dfield_add_member(struct cow_dfield *f, const char *name)
{
  if (f->committed) return;
  f->n_members++;
  f->members = (char**) realloc(f->members, f->n_members*sizeof(char*));
  f->members[f->n_members-1] = (char*) malloc(strlen(name)+1);
  strcpy(f->members[f->n_members-1], name);
}
const char *cow_dfield_iteratemembers(struct cow_dfield *f)
{
  f->member_iter = 0;
  return cow_dfield_nextmember(f);
}
const char *cow_dfield_nextmember(struct cow_dfield *f)
{
  return f->member_iter++ < f->n_members ? f->members[f->member_iter-1] : NULL;
}
void cow_dfield_commit(struct cow_dfield *f)
{
  if (f->committed) return;
  const int n_zones = cow_domain_getnumlocalzones(f->domain);
  f->data = malloc(n_zones * f->n_members * sizeof(double));
  f->committed = 1;
}

int main()
{
  struct cow_domain *domain = cow_domain_new();
  struct cow_dfield *prim = cow_domain_addfield(domain, "primitive");
  struct cow_dfield *magf = cow_domain_addfield(domain, "magnetic");

#if (COW_MPI)
  printf("was compiled with COW_MPI\n");
#endif

  cow_dfield_add_member(prim, "vx");
  cow_dfield_add_member(prim, "vy");
  cow_dfield_add_member(prim, "vz");

  cow_dfield_add_member(magf, "Bx");
  cow_dfield_add_member(magf, "By");
  cow_dfield_add_member(magf, "Bz");

  cow_domain_commit(domain);

  for (struct cow_dfield *d = cow_domain_iteratefields(domain);
       d != NULL; d = cow_domain_nextfield(domain)) {
    printf("%s\n", cow_dfield_getname(d));
    for (const char *m = cow_dfield_iteratemembers(d);
	 m != NULL; m = cow_dfield_nextmember(d)) {
      printf("\t%s\n", m);
    }
  }

  cow_domain_del(domain);
  return 0;
}
