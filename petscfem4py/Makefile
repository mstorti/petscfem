.PHONY: default srd build install clean

default: build

src:
	python setup.py build_src --force
build:
	python setup.py build

install:
	python setup.py install --home=${HOME}

clean:
	python setup.py clean --all
