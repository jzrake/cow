

.PHONY : clib python cython

default : clib

all : python clib

clib : 
	@make -C src

python :
	@make -C python

cython : clib
	@make -C cython

clean :
	@make -C src clean
	@make -C python clean
	@make -C cython clean
	@rm -rf lib bin include
