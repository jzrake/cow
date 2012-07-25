
COWPY_INSTALL ?= $(PWD)

.PHONY: COWPY

default : COWPY

COWPY : 
	python setup.py build
	python setup.py install --prefix=$(COWPY_INSTALL)

clean :
	python setup.py clean --all

