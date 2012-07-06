
CC           ?= mpicc
CXX          ?= mpicxx
COW_HDF5     ?= 0
COW_MPI      ?= 0
COW_HDF5_MPI ?= 0
COW_FFTW     ?= 0
CFLAGS       ?= -Wall -g -O0
HDF5_HOME    ?= /usr/local
FFTW_HOME    ?= /usr/local

ifeq ($(COW_HDF5), 1)
INC += -I$(HDF5_HOME)/include
LIB += -L$(HDF5_HOME)/lib -lz -lhdf5
endif

ifeq ($(COW_FFTW), 1)
INC += -I$(FFTW_HOME)/include
LIB += -L$(FFTW_HOME)/lib -lfftw
endif

DEFINES = \
	-DCOW_MPI=$(COW_MPI) \
	-DCOW_HDF5=$(COW_HDF5) \
	-DCOW_HDF5_MPI=$(COW_HDF5_MPI) \
	-DCOW_FFTW=$(COW_FFTW) \
	-DFFT_FFTW

OBJ = cow.o hist.o io.o fft.o fft_3d.o pack_3d.o remap_3d.o factor.o
EXE = main makehist testhist milos mhdstats

default : $(EXE)

%.o : %.c
	$(CC) $(CFLAGS) -o $@ $< $(DEFINES) $(INC) -c -std=c99

%.o : %.cpp
	$(CC) $(CFLAGS) -o $@ $< $(DEFINES) $(INC) -c

main : main.o $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LIB)

testhist : testhist.o $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LIB)

makehist : makehist.o $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LIB)

mhdstats : mhdstats.o $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LIB)

milos : milos.o $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LIB)

clean :
	rm -rf $(EXE) $(OBJ)
