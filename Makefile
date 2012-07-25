

# install prefix (use as command line variable)
PREFIX ?= /usr/local

# install prefix (use as environment variable)
COWLIB_INSTALL ?= $(PREFIX)

.PHONY : clib python cython install

default : clib

all : python clib

clib : 
	@make -C src

python :
	@make -C python

cython : clib
	@make -C cython

install : clib
	mkdir -p $(COWLIB_INSTALL)/include; cp include/* $(COWLIB_INSTALL)/include
	mkdir -p $(COWLIB_INSTALL)/bin; cp bin/* $(COWLIB_INSTALL)/bin
	mkdir -p $(COWLIB_INSTALL)/lib; cp lib/* $(COWLIB_INSTALL)/lib

clean :
	@make -C src clean
	@make -C python clean
	@make -C cython clean
	@rm -rf lib bin include
