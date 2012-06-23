

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


struct cow_dfield; // forward declarations (for opaque data structure)
struct cow_domain;

struct cow_dfield *cow_domain_add_field(struct cow_domain *d, const char *name);
struct cow_dfield *cow_domain_get_field(struct cow_domain *d, const char *name);
struct cow_dfield *cow_dfield_new();
struct cow_domain *cow_domain_new();
void cow_dfield_del(struct cow_dfield *f);
void cow_dfield_add_member(struct cow_dfield *f, const char *name);
void cow_dfield_set_name(struct cow_dfield *f, const char *name);
const char *cow_dfield_get_name(struct cow_dfield *f);

void cow_domain_set_size(struct cow_domain *d, int dim, int size);
void cow_domain_set_ndim(struct cow_domain *d, int ndim);
void cow_domain_set_guard(struct cow_domain *d, int guard);


// -----------------------------------------------------------------------------
//
// struct cow_domain interface functions
//
// -----------------------------------------------------------------------------
struct cow_domain
{
  int A_nint[3];
  int L_ntot[3];
  int L_strt[3];
  int G_ntot[3];
  int G_strt[3];
  int n_dims;
  int n_ghst;
  int n_fields;
  int committed;
  struct cow_dfield **fields;
} ;

struct cow_domain *cow_domain_new()
{
  struct cow_domain dom = {
    .A_nint = { 1, 1, 1 },
    .L_ntot = { 1, 1, 1 },
    .L_strt = { 1, 1, 1 },
    .G_ntot = { 1, 1, 1 },
    .G_strt = { 1, 1, 1 },
    .n_dims = 1,
    .n_ghst = 1,
    .n_fields = 0,
    .committed = 0,
    .fields = NULL
  } ;
  struct cow_domain *d = (struct cow_domain*) malloc(sizeof(struct cow_domain));
  *d = dom;
  return d;
}
void cow_domain_del(struct cow_domain *d)
{
  for (int n=0; n<d->n_fields; ++n) cow_dfield_del(d->fields[n]);
  free(d->fields);
  free(d);
}
struct cow_dfield *cow_domain_add_field(struct cow_domain *d, const char *name)
{
  d->n_fields++;
  d->fields = (struct cow_dfield**) realloc(d->fields,
					    d->n_fields*sizeof(struct cow_dfield*));
  struct cow_dfield *f = cow_dfield_new();
  cow_dfield_set_name(f, name);
  return d->fields[d->n_fields-1] = f;
}
struct cow_dfield *cow_domain_get_field(struct cow_domain *d, const char *name)
{
  for (int n=0; n<d->n_fields; ++n) {
    if (strcmp(cow_dfield_get_name(d->fields[n]), name) == 0) {
      return d->fields[n];
    }
  }
  return NULL;
}
void cow_domain_set_size(struct cow_domain *d, int dim, int size)
{
  if (dim > 3) return;
  d->G_ntot[dim] = size;
}
void cow_domain_set_ndim(struct cow_domain *d, int ndim)
{
  if (ndim > 3) return;
  d->n_dims = ndim;
}
void cow_domain_set_guard(struct cow_domain *d, int guard)
{
  d->n_ghst = guard;
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
  int n_members;
} ;
struct cow_dfield *cow_dfield_new()
{
  struct cow_dfield *f = (struct cow_dfield*) malloc(sizeof(struct cow_dfield));
  f->name = NULL;
  f->members = NULL;
  f->n_members = 0;
  return f;
}
void cow_dfield_del(struct cow_dfield *f)
{
  for (int n=0; n<f->n_members; ++n) free(f->members[n]);
  free(f->members);
  free(f->name);
  free(f);
}
void cow_dfield_set_name(struct cow_dfield *f, const char *name)
{
  f->name = (char*) realloc(f->name, strlen(name)+1);
  strcpy(f->name, name);
}
const char *cow_dfield_get_name(struct cow_dfield *f)
{
  return f->name;
}
void cow_dfield_add_member(struct cow_dfield *f, const char *name)
{
  f->n_members++;
  f->members = (char**) realloc(f->members, f->n_members*sizeof(char*));
  f->members[f->n_members-1] = (char*) malloc(strlen(name)+1);
  strcpy(f->members[f->n_members-1], name);
}



int main()
{
  struct cow_domain *domain = cow_domain_new();
  struct cow_dfield *prim = cow_domain_add_field(domain, "primitive");

  cow_dfield_add_member(prim, "vx");
  cow_dfield_add_member(prim, "vy");

  for (int n=0; n<prim->n_members; ++n) {
    printf("%d: %s\n", n, prim->members[n]);
  }

  cow_domain_del(domain);
  return 0;
}
