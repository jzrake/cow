
print "Hello world!"

cimport numpy as np

cdef extern from "cow.h":
    cdef enum:
        COW_NOREOPEN_STDOUT = (2<<0)
        COW_NOREOPEN_STDOUT = (2<<0)
        COW_DISABLE_MPI     = (2<<1)
        COW_ALL_DIMS             = -41
        COW_HIST_SPACING_LINEAR  = -42
        COW_HIST_SPACING_LOG     = -43
        COW_HIST_BINMODE_COUNTS  = -44 # traditional histogram
        COW_HIST_BINMODE_DENSITY = -45 # divides by bin width
        COW_HIST_BINMODE_AVERAGE = -46
        COW_PROJECT_OUT_DIV      = -47 # used for Helmholtz decomposition
        COW_PROJECT_OUT_CURL     = -48
        COW_SAMPLE_NEAREST       = -49 # sample the nearest zone center
        COW_SAMPLE_LINEAR        = -50 # use (uni/bi/tri) linear interp
        COW_SAMPLE_ERROR_OUT     = -51 # out-of-bounds sample request
        COW_SAMPLE_ERROR_WRONGD  = -52 # wrong number of dims on sample coords

    cdef struct cow_domain
    cdef struct cow_dfield
    cdef struct cow_histogram
    ctypedef void (*cow_transform)(double *result, double **args, int **strides,
                                   void *udata)

    void cow_init(int argc, char **argv, int modes)
    void cow_finalize()
    int cow_mpirunning()

    cow_domain *cow_domain_new()
    void cow_domain_commit(cow_domain *d)
    void cow_domain_del(cow_domain *d)
    void cow_domain_setsize(cow_domain *d, int dim, int size)
    void cow_domain_setndim(cow_domain *d, int ndim)
    void cow_domain_setguard(cow_domain *d, int guard)
    void cow_domain_setprocsizes(cow_domain *d, int dim, int size)
    void cow_domain_setcollective(cow_domain *d, int mode)
    void cow_domain_setchunk(cow_domain *d, int mode)
    void cow_domain_setalign(cow_domain *d, int alignthreshold, int diskblocksize)
    void cow_domain_readsize(cow_domain *d, char *fname, char *dname)
    int cow_domain_getndim(cow_domain *d)
    int cow_domain_getguard(cow_domain *d)
    int cow_domain_getnumlocalzonesincguard(cow_domain *d, int dim)
    int cow_domain_getnumlocalzonesinterior(cow_domain *d, int dim)
    int cow_domain_getnumglobalzones(cow_domain *d, int dim)
    int cow_domain_getglobalstartindex(cow_domain *d, int dim)
    int cow_domain_getgridspacing(cow_domain *d, int dim)
    int cow_domain_getcartrank(cow_domain *d)
    int cow_domain_getcartsize(cow_domain *d)
    int cow_domain_subgridatposition(cow_domain *d, double x, double y, double z)
    int cow_domain_indexatposition(cow_domain *d, int dim, double x)
    double cow_domain_positionatindex(cow_domain *d, int dim, int index)
    void cow_domain_barrier(cow_domain *d)

    cow_dfield *cow_dfield_new(cow_domain *domain, char *name)
    cow_dfield *cow_dfield_dup(cow_dfield *f)
    void cow_dfield_commit(cow_dfield *f)
    void cow_dfield_del(cow_dfield *f)
    void cow_dfield_addmember(cow_dfield *f, char *name)
    void cow_dfield_setname(cow_dfield *f, char *name)
    void cow_dfield_extract(cow_dfield *f, int *I0, int *I1, void *out)
    void cow_dfield_replace(cow_dfield *f, int *I0, int *I1, void *out)
    void cow_dfield_loop(cow_dfield *f, cow_transform op, void *udata)
    void cow_dfield_settransform(cow_dfield *f, cow_transform op)
    void cow_dfield_clearargs(cow_dfield *f)
    void cow_dfield_pusharg(cow_dfield *f, cow_dfield *arg)
    void cow_dfield_setuserdata(cow_dfield *f, void *userdata)
    void cow_dfield_setiparam(cow_dfield *f, int p)
    void cow_dfield_setfparam(cow_dfield *f, double p)
    void cow_dfield_transformexecute(cow_dfield *f)
    char *cow_dfield_iteratemembers(cow_dfield *f)
    char *cow_dfield_nextmember(cow_dfield *f)
    char *cow_dfield_getname(cow_dfield *f)
    cow_domain *cow_dfield_getdomain(cow_dfield *f)
    int cow_dfield_getstride(cow_dfield *f, int dim)
    int cow_dfield_getnmembers(cow_dfield *f)
    size_t cow_dfield_getdatabytes(cow_dfield *f)
    void cow_dfield_setbuffer(cow_dfield *f, void *buffer)
    void cow_dfield_sampleglobalind(cow_dfield *f, int i, int j, int k, double **x, int *n0)
    int cow_dfield_setsamplecoords(cow_dfield *f, double *x, int n0, int n1)
    void cow_dfield_getsamplecoords(cow_dfield *f, double **x, int *n0, int *n1)
    void cow_dfield_getsampleresult(cow_dfield *f, double **x, int *n0, int *n1)
    void cow_dfield_setsamplemode(cow_dfield *f, int mode)
    void cow_dfield_sampleexecute(cow_dfield *f)
    int cow_dfield_getownsdata(cow_dfield *f)
    void *cow_dfield_getbuffer(cow_dfield *f)
    void cow_dfield_syncguard(cow_dfield *f)
    void cow_dfield_reduce(cow_dfield *f, double x[3])
    void cow_dfield_write(cow_dfield *f, char *fname)
    void cow_dfield_read(cow_dfield *f, char *fname)

    cow_histogram *cow_histogram_new()
    void cow_histogram_commit(cow_histogram *h)
    void cow_histogram_del(cow_histogram *h)
    void cow_histogram_setbinmode(cow_histogram *h, int binmode)
    void cow_histogram_setspacing(cow_histogram *h, int spacing)
    void cow_histogram_setnbins(cow_histogram *h, int dim, int nbinsx)
    void cow_histogram_setlower(cow_histogram *h, int dim, double v0)
    void cow_histogram_setupper(cow_histogram *h, int dim, double v1)
    void cow_histogram_setfullname(cow_histogram *h, char *fullname)
    void cow_histogram_setnickname(cow_histogram *h, char *nickname)
    void cow_histogram_setdomaincomm(cow_histogram *h, cow_domain *d)
    void cow_histogram_addsample1(cow_histogram *h, double x, double w)
    void cow_histogram_addsample2(cow_histogram *h, double x, double y, double w)
    void cow_histogram_dumpascii(cow_histogram *h, char *fn)
    void cow_histogram_dumphdf5(cow_histogram *h, char *fn, char *dn)
    void cow_histogram_seal(cow_histogram *h)
    int cow_histogram_getsealed(cow_histogram *h)
    long cow_histogram_gettotalcounts(cow_histogram *h)
    void cow_histogram_populate(cow_histogram *h, cow_dfield *f, cow_transform op)
    void cow_histogram_getbinlocx(cow_histogram *h, double **x, int *n0)
    void cow_histogram_getbinlocy(cow_histogram *h, double **x, int *n0)
    void cow_histogram_getbinval1(cow_histogram *h, double **x, int *n0)
    void cow_histogram_getbinval2(cow_histogram *h, double **x, int *n0, int *n1)
    double cow_histogram_getbinval(cow_histogram *h, int i, int j)
    char *cow_histogram_getname(cow_histogram *h)

    void cow_fft_pspecvecfield(cow_dfield *f, cow_histogram *h)
    void cow_fft_helmholtzdecomp(cow_dfield *f, int mode)



cdef class DistrbutedDomain(object):
    cdef cow_domain *_c
    def __cinit__(self):
        print "building domain"
        self._c = cow_domain_new()

    def __dealloc__(self):
        print "killing domain"
        cow_domain_del(self._c)
