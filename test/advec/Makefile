#__INSERT_LICENSE__
# $Id: Makefile,v 1.8 2002/12/15 18:23:26 mstorti Exp $ 

.PHONY: force all run default

default: rotating_cone
tests: skew parallel_cone skew_cone rotating_cone

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include ../Makefile.base.test
SRCS = 

#s Linear Advection tests
LOC_PROG := advdif
CASE_NAME := advec
PROG := $(ADVDIF_O)

skew:
	$(MAKE) case_name=$@ run
	octave -qH procskew.m > $@.verif.tmp

cone:
	$(MAKE) case_name=$@ run
	octave -qH proc.m > $@.verif.tmp

parallel_cone:
	$(MAKE) case_name=$@ run
	octave -qH verif_$@.m > $@.verif.tmp

skew_cone:
	$(MAKE) case_name=$@ run
	octave -qH verif_$@.m > $@.verif.tmp

rotating_cone:
	$(MAKE) case_name=$@ run
	octave -qH verif_$@.m > $@.verif.tmp

advec.depl: force

case := $(CASE_NAME)
NP := $(shell procsel $(PETSCFEM_PROCTABLE))
run: advec.depl
	-$(MPIRUN) -np 1 -machinefile machi.dat $(PROG) \
		-case $(case).depl  > advec.$(case_name).out.tmp

dx: viewero.dx
viewero.dx: viewer.dx
	eperl -P $< > $@

local_clean::

#s
