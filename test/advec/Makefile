#__INSERT_LICENSE__
# $Id: Makefile,v 1.1 2002/12/13 16:13:54 mstorti Exp $ 

.PHONY: force all run 

all: run

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include ../Makefile.base.test
SRCS = 

#s Linear Advection tests
LOC_PROG := advdif
CASE_NAME := advec
PROG := $(ADVDIF)

advec.depl: force

case := $(CASE_NAME)
NP := $(shell procsel $(PETSCFEM_PROCTABLE))
run: advec.depl
	-$(MPIRUN) -np 1 -machinefile machi.dat $(ADVDIF) \
			-case $(case).depl # > $(case).out

local_clean::

#s
