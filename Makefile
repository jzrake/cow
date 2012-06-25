
CC = mpicc
COW_HDF5 = 1
COW_MPI = 1
COW_HDF5_MPI = 1
CFLAGS = -Wall -g -O0 -std=c99

HDF5_HOME ?= usr/local
INC = -I$(HDF5_HOME)/include
LIB = -L$(HDF5_HOME)/lib -lz -lhdf5

DEFINES = \
	-DCOW_MPI=$(COW_MPI) \
	-DCOW_HDF5=$(COW_HDF5) \
	-DCOW_HDF5_MPI=$(COW_HDF5_MPI)

default : main

%.o : %.c
	$(CC) $(CFLAGS) -o $@ $< $(DEFINES) $(INC) -c

main : main.o cow.o h5mpi.o
	$(CC) $(CFLAGS) -o $@ $^ $(DEFINES) $(LIB)

clean :
	rm -rf main *.o *.dSYM
