# $Id: Makefile,v 1.1.2.4 2006/05/30 20:38:33 dalcinl Exp $

.PHONY: default src build install uninstall clean swig

default: build

src:
	python setup.py build_src

build:
	python setup.py build

install:
	python setup.py install --home=${HOME}

uninstall:
	${RM} -r ${HOME}/lib/python/petscfem4py

clean:
	python setup.py clean --all

swig:
	python setup.py build_src --force

