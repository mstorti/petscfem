#__INSERT_LICENSE__
# $Id: Makefile,v 1.5 2002/12/14 21:19:05 mstorti Exp $ 

.PHONY: force all run default

default: skew

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include ../Makefile.base.test
SRCS = 

#s Linear Advection tests
LOC_PROG := advdif
CASE_NAME := advec
PROG := $(ADVDIF_G)

skew:
	$(MAKE) case_name=skew run

cone:
	$(MAKE) case_name=cone run

rotating_cone:
	$(MAKE) case_name=rotating_cone run

advec.depl: force

case := $(CASE_NAME)
NP := $(shell procsel $(PETSCFEM_PROCTABLE))
run: advec.depl
	-$(MPIRUN) -np 1 -machinefile machi.dat $(PROG) \
		-case $(case).depl # > $(case).out

dx: viewero.dx
viewero.dx: viewer.dx
	eperl -P $< > $@

local_clean::

#s
