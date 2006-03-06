# $Id: Makefile,v 1.1.2.3 2006/03/06 16:56:04 rodrigop Exp $

.PHONY: default src swig  build install clean

default: build

src:
	python setup.py build_src

swig:
	python setup.py build_src --force

build:
	python setup.py build

install:
	python setup.py install --home=${HOME}

clean:
	python setup.py clean --all
