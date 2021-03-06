
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
SRC = cow.c fft.c fft_3d.c hist.c \
	io.c pack_3d.c remap_3d.c samp.c \
	read-ff-sdf.c

OBJ = $(SRC:.c=.o)
DEP = $(SRC:.c=.dep)
LIB = libcow.a
EXE = ffe-pspec2d


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
default : $(EXE)

$(EXE) : $(EXE).c $(LIB)
	$(CC) $(CFLAGS) $^ $(HDF5_I) $(FFTW_I) $(HDF5_L) $(FFTW_L) $(RNPL_L) $(CLIBS) -o $@

$(LIB) : $(OBJ)
	$(AR) $@ $?
	$(RANLIB) $@

%.o : %.c cow-cfg.h $(MAKEFILE_IN)
	$(CCMM) -MM $< > $(<:.c=.dep) $(HDF5_I) $(FFTW_I) $(RNPL_I)
	$(CC) $(CFLAGS) -c $< $(HDF5_I) $(FFTW_I) $(RNPL_I)

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
	@echo "OBJ: $(OBJ)"
	@echo "SRC: $(SRC)"
	@echo "DEP: $(DEP)"

clean :
	$(RM) $(LIB) $(OBJ) $(DEP) $(EXE)

.FORCE :
.PHONY : show clean

-include *.dep
