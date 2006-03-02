.PHONY:  default all seq mpi

export PYTHONPATH=${HOME}/lib/python

default: seq

all: seq mpi

PYTHON=python

seq:
	@for test_py in `ls test_*.py`; do \
          echo runing $$test_py ... && ${PYTHON} $$test_py; \
	done

MPIPYTHON=python
MPIEXEC=mpiexec -n
NP=2

mpi:
	@for test_py in `ls test_*.py`; do \
	  ${MPIEXEC} ${NP} ${MPIPYTHON} $$test_py; \
	done
