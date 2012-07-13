
CC           ?= mpicc
COW_HDF5     ?= 0
COW_MPI      ?= 0
COW_HDF5_MPI ?= 0
COW_FFTW     ?= 0
CFLAGS       ?= -Wall -g -O0
FFTW_INC     ?= 
FFTW_LIB     ?= 
HDF5_INC     ?= 
HDF5_LIB     ?= 


DEFINES = \
	-DCOW_MPI=$(COW_MPI) \
	-DCOW_HDF5=$(COW_HDF5) \
	-DCOW_HDF5_MPI=$(COW_HDF5_MPI) \
	-DCOW_FFTW=$(COW_FFTW)

LIB = $(HDF5_LIB) $(FFTW_LIB)
INC = $(HDF5_INC) $(FFTW_INC)

OBJ = cow.o hist.o io.o samp.o fft.o fft_3d.o pack_3d.o remap_3d.o factor.o
EXE = main makehist milos mhdstats srhdhist testhist testfft testsamp

.PHONY: COWPY

default : COWPY

COWPY : 
	python setup.py build
	python setup.py install --prefix=$(COW_INSTALL)

%.o : %.c
	$(CC) $(CFLAGS) -o $@ $< $(DEFINES) $(INC) -c -std=c99

cow_wrap.c : cow.i
	swig -python -o $@ $^

cowpy.o : cowpy.c
	$(CC) $(CFLAGS) -o $@ $< -c -I$(PY_INC) $(DEFINES) -std=c99

cow_wrap.o : cow_wrap.c
	$(CC) $(CFLAGS) -o $@ $< -c -I$(PY_INC)

cowpy : cowpy.o cow_wrap.o $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LIB) -l$(PY_LIB)

main : main.o $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LIB)

makehist : makehist.o $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LIB)

mhdstats : mhdstats.o $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LIB)

milos : milos.o $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LIB)

srhdhist : srhdhist.o $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LIB)

testhist : testhist.o $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LIB)

testfft : testfft.o $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LIB)

testsamp : testsamp.o $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LIB)

clean :
	rm -rf $(EXE) $(OBJ) cow_wrap.* cow.py lib
	python setup.py clean --all
