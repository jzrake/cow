
CC ?= mpicc
CXX ?= mpicxx
COW_HDF5 ?= 1
COW_MPI ?= 1
COW_HDF5_MPI ?= 1
CFLAGS ?= -Wall -g -O0

HDF5_HOME ?= /usr/local
INC = -I$(HDF5_HOME)/include
LIB = -L$(HDF5_HOME)/lib -lz -lhdf5

DEFINES = \
	-DCOW_MPI=$(COW_MPI) \
	-DCOW_HDF5=$(COW_HDF5) \
	-DCOW_HDF5_MPI=$(COW_HDF5_MPI)

OBJ = cow.o hist.o fft_3d.o io.o pack_3d.o remap_3d.o

EXE = main makehist testhist milos cowpy

default : $(EXE)


%.o : %.c
	$(CC) $(CFLAGS) -o $@ $< $(DEFINES) $(INC) -c -std=c99

cow_wrap.c : cow.i
	swig -python -o $@ $^

cowpy.o : cowpy.c
	$(CC) $(CFLAGS) -o $@ $< -c $(DEFINES) -I/Library/Frameworks/Python.framework/Headers -std=c99

cow_wrap.o : cow_wrap.c
	$(CC) $(CFLAGS) -o $@ $< -c -I/Library/Frameworks/Python.framework/Headers

cowpy : cowpy.o cow_wrap.o $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(DEFINES) $(LIB) -lpython

main : main.o $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(DEFINES) $(LIB)

testhist : testhist.o $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(DEFINES) $(LIB)

makehist : makehist.o $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(DEFINES) $(LIB)

milos : milos.o $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(DEFINES) $(LIB)

clean :
	rm -rf $(EXE) $(OBJ) cow_wrap.* cow.py

