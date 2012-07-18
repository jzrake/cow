
COW_INSTALL ?= $(PWD)

.PHONY: COWPY

default : COWPY

COWPY : 
	python setup.py build
	python setup.py install --prefix=$(COW_INSTALL)

clean :
	python setup.py clean --all
	rm -rf cow.py cow_wrap.c *.pyc lib
