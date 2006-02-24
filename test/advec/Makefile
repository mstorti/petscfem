#__INSERT_LICENSE__
# $Id: Makefile,v 1.21.46.2 2006/02/24 12:26:57 rodrigop Exp $ 

.PHONY: force all run default skew parallel_cone skew_cone rotating_cone	\
	skew_noise parallel_cone_noise skew_cone_noise rotating_cone_noise	\
	runs

default: parallel_cone
tests: tests_noise0 tests_noise1

tests_noise0:  skew parallel_cone skew_cone rotating_cone

tests_noise1: skew_noise parallel_cone_noise skew_cone_noise rotating_cone_noise

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include ../Makefile.base.test
SRCS = 

#s Linear Advection tests
LOC_PROG := advdif
# CASE_NAME := smoke
CASE_NAME := advec
PROG := $(ADVDIF)
# PROG := $(NS)

skew:
	$(MAKE) case_name=$@ run
	$(OCTAVE) verif_$@.m > $@.verif.tmp
#	$(OCTAVE) verif_$@.m

skew_noise:
	$(MAKE) case_name=$@ run
	$(OCTAVE) verif_$@.m > $@.verif.tmp

# cone:
# 	$(MAKE) case_name=$@ run
# 	$(OCTAVE) proc.m > $@.verif.tmp

parallel_cone:
	$(MAKE) case_name=$@ run
	$(OCTAVE) verif_$@.m > $@.verif.tmp

parallel_cone_noise:
	$(MAKE) case_name=$@ run
	$(OCTAVE) verif_$@.m > $@.verif.tmp

skew_cone:
	$(MAKE) case_name=$@ run
	$(OCTAVE) verif_$@.m > $@.verif.tmp

skew_cone_noise:
	$(MAKE) case_name=$@ run
	$(OCTAVE) verif_$@.m > $@.verif.tmp

rotating_cone:
	$(MAKE) case_name=$@ run
	$(OCTAVE) verif_$@.m > $@.verif.tmp

rotating_cone_noise:
	$(MAKE) case_name=$@ run
	$(OCTAVE) verif_$@.m > $@.verif.tmp

advec.depl: force

case := $(CASE_NAME)
NP := $(shell procsel $(PETSCFEM_PROCTABLE))
run: advec.depl
	-$(MPIRUN) -np 2     $(PROG) \
		-case $(case).depl  > advec.$(case_name).out.tmp

runs: CASE_NAME := smoke
runs: smoke.depl
	-$(MPIRUN) -np 1     $(PROG) \
		-case $(CASE_NAME).depl >  $(case_name).verif.tmp

smoke.depl: force

dx: viewero.dx
viewero.dx: viewer.dx
	$(EPERL) -P $< > $@

local_clean::

#s
