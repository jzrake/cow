
# if there is no Makefile.in then use the template
# --------------------------------------------------
ifeq ($(shell test -e Makefile.in && echo 1), 1)
MAKEFILE_IN = Makefile.in
else
MAKEFILE_IN =  Makefile.in.template
endif
include $(MAKEFILE_IN)


GIT_SHA := $(shell git rev-parse HEAD | cut -c 1-10)
TMP_H := $(shell mktemp -u make.XXXXXX)


# object code required for executables
# --------------------------------------------------
SRC = $(wildcard *.c)
OBJ = $(SRC:.c=.o)
DEP = $(OBJ:.o=.dep)
LIB = libcow.a
EXE = ffe-pspec2d ffe


ifeq ($(HAVE_HDF5), 1)
HDF5_I = -I$(HDF5_HOME)/include
HDF5_L = -L$(HDF5_HOME)/lib -lhdf5
endif

ifeq ($(HAVE_FFTW), 1)
FFTW_I = -I$(FFTW_HOME)/include
FFTW_L = -L$(FFTW_HOME)/lib -lfftw3
endif

ifeq ($(HAVE_RNPL), 1)
RNPL_I = -I$(RNPL_HOME)/include
RNPL_L = -L$(RNPL_HOME)/lib -lrnpl
endif


# build rules
# --------------------------------------------------
default : ffe ffe-pspec2d

ffe-pspec2d : ffe-pspec2d.o $(LIB)
	$(CC) $(CFLAGS) $(CLIBS) $^ $(HDF5_L) $(FFTW_L) $(RNPL_L) -o $@

ffe : ffe.o $(LIB)
	$(CC) $(CFLAGS) $(CLIBS) $^ $(HDF5_L) $(FFTW_L) $(RNPL_L) -o $@

$(LIB) : $(OBJ)
	$(AR) $@ $?
	$(RANLIB) $@

%.o : %.c cow-cfg.h $(MAKEFILE_IN)
	@$(CC) -MM $< > $(<:.c=.dep) $(HDF5_I) $(FFTW_I) $(RNPL_I)
	$(CC) $(CFLAGS) $< -c $(HDF5_I) $(FFTW_I) $(RNPL_I)

cow-cfg.h : .FORCE
	@echo "/* cow config header file */" > $(TMP_H)
	@echo "#define COW_MPI $(HAVE_MPI)" >> $(TMP_H)
	@echo "#define COW_HDF5 $(HAVE_HDF5)" >> $(TMP_H)
	@echo "#define COW_FFTW $(HAVE_FFTW)" >> $(TMP_H)
	@echo "#define COW_RNPL $(HAVE_RNPL)" >> $(TMP_H)
	@echo "#define COW_GIT_SHA \"$(GIT_SHA)\"" >> $(TMP_H)
	@cmp -s $(TMP_H) $@ || (echo "[cow-cfg.h updated]"; cat $(TMP_H)>$@)
	@$(RM) $(TMP_H)

show :
	@echo "MAKEFILE_IN: $(MAKEFILE_IN)"
	@echo "CC: $(CC)"
	@echo "CFLAGS: $(CFLAGS)"
	@echo "EXE: $(EXE)"
	@echo "SRC: $(SRC)"
	@echo "OBJ: $(OBJ)"
	@echo "DEP: $(DEP)"

clean :
	$(RM) $(LIB) $(OBJ) $(DEP) $(EXE)

.FORCE :
.PHONY : show clean

-include *.dep
