
COW_INSTALL ?= $(PWD)

.PHONY: COWPY

default : COWPY

COWPY : 
	python setup.py build
	python setup.py install --prefix=$(COW_INSTALL)

clean :
	python setup.py clean --all

