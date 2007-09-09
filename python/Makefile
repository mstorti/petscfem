# $Id$

.PHONY: default src build install uninstall clean swig

default: build

src:
	python setup.py build_src

build:
	python setup.py build

install:
	python setup.py install --home=${HOME}

uninstall:
	${RM} -r ${HOME}/lib/python/pf4py

clean:
	python setup.py clean --all

swig:
	python setup.py build_src --force


rebuild:
	rm -rf build/lib*/pf4py/kernel/_core.so
	make -C .. ns_all
	make
