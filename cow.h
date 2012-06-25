

#ifndef COW_HEADER_INCLUDED
#define COW_HEADER_INCLUDED
#include <stdlib.h>
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
void cow_domain_setsize(cow_domain *d, int dim, int size);
int cow_domain_getsize(cow_domain *d, int dim);
void cow_domain_setndim(cow_domain *d, int ndim);
void cow_domain_setguard(cow_domain *d, int guard);
int cow_domain_getguard(cow_domain *d);
void cow_domain_setprocsizes(cow_domain *d, int dim, int size);
//cow_dfield *cow_domain_addfield(cow_domain *d, const char *name);
//cow_dfield *cow_domain_iteratefields(cow_domain *d);
//cow_dfield *cow_domain_nextfield(cow_domain *d);
int cow_domain_getnumlocalzones(cow_domain *d);

cow_dfield *cow_dfield_new(cow_domain *domain, const char *name);
void cow_dfield_commit(cow_dfield *f);
void cow_dfield_del(cow_dfield *f);
void cow_dfield_addmember(cow_dfield *f, const char *name);
void cow_dfield_setname(cow_dfield *f, const char *name);
void cow_dfield_extract(cow_dfield *f, const int *I0, const int *I1, void *out);
void cow_dfield_replace(cow_dfield *f, const int *I0, const int *I1, void *out);

int cow_dfield_getstride(cow_dfield *d, int dim);
const char *cow_dfield_getname(cow_dfield *f);
const char *cow_dfield_iteratemembers(cow_dfield *f);
const char *cow_dfield_nextmember(cow_dfield *f);
void *cow_dfield_getdata(cow_dfield *f);
void cow_dfield_syncguard(cow_dfield *f);

#endif // COW_HEADER_INCLUDED
