
.PHONY : python

default : clib

clib :
	make -C src

python :
	make -C python

all : python clib

clean :
	make -C src clean
	make -C python clean
