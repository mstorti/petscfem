#__INSERT_LICENSE__
# $Id: Makefile,v 1.3 2002/12/14 14:04:38 mstorti Exp $ 

.PHONY: force all run default

default: cone

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include ../Makefile.base.test
SRCS = 

#s Linear Advection tests
LOC_PROG := advdif
CASE_NAME := advec
PROG := $(ADVDIF_O)

cone:
	$(MAKE) case_name=cone run

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
