#__INSERT_LICENSE__

.PHONY: force all run default skew parallel_cone skew_cone rotating_cone	\
	skew_noise parallel_cone_noise skew_cone_noise rotating_cone_noise	\
	runs

default: parallel_cone
tests: tests_noise0 tests_noise1

tests_noise0:  skew parallel_cone skew_cone rotating_cone

tests_noise1: skew_noise parallel_cone_noise skew_cone_noise \
		 rotating_cone_noise

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
# MCH :=  -machinefile machi.dat 
run: advec.depl
	-$(MPIEXEC) $(MCH) -n 2 $(PROG) \
		-case $(case).depl  $(call redir,advec.$(case_name).out.tmp)

runs: CASE_NAME := smoke
runs: smoke.depl
	-$(MPIEXEC) $(MCH) -n 1 $(PROG) \
		-case $(CASE_NAME).depl $(call in_backg,$@)

smoke.depl: force

dx: viewero.dx
viewero.dx: viewer.dx
	$(EPERL) -P $< > $@

local_clean::

prog:
	$(MAKE) -C $(PETSCFEM_DIR) advdif

#s

