#__INSERT_LICENSE__
# $Id: Makefile,v 1.6 2002/12/15 04:01:09 mstorti Exp $ 

.PHONY: force all run default

default: rotating_cone

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include ../Makefile.base.test
SRCS = 

#s Linear Advection tests
LOC_PROG := advdif
CASE_NAME := advec
PROG := $(ADVDIF_O)

skew:
	$(MAKE) case_name=skew run

cone:
	$(MAKE) case_name=cone run

parallel_cone:
	$(MAKE) case_name=parallel_cone run

rotating_cone:
	$(MAKE) case_name=rotating_cone run

skew_cone:
	$(MAKE) case_name=skew_cone run

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
