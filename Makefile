

# install prefix (use as command line variable)
PREFIX ?= $(pwd)

# install prefix (use as environment variable)
COWLIB_INSTALL ?= $(PREFIX)

.PHONY : clib install cowpy

default : clib

all : clib cowpy

clib : 
	@make -C src

cowpy : clib
	@make -C cowpy

install : clib
	mkdir -p $(COWLIB_INSTALL)/include; cp include/* $(COWLIB_INSTALL)/include
	mkdir -p $(COWLIB_INSTALL)/bin; cp bin/* $(COWLIB_INSTALL)/bin
	mkdir -p $(COWLIB_INSTALL)/lib; cp lib/* $(COWLIB_INSTALL)/lib

clean :
	@make -C src clean
	@make -C cowpy clean
	@rm -rf lib bin include
