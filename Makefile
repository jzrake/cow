

.PHONY : clib python

default : clib

all : python clib

clib : 
	@make -C src

python :
	@make -C python

clean :
	@make -C src clean
#	@make -C python clean
	@rm -rf lib bin include
