#__INSERT_LICENSE__
# $Id: Makefile,v 1.9 2003/02/13 21:57:52 mstorti Exp $ 

.PHONY: force all run tests data scalar diag full constraint default

default: constraint

tests: scalar diag full constraint fstack

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test
SRCS = 

LOC_PROG := ns
CASE_NAME := genload
PROG := $(NS)

oct_check := yes
export oct_check

.PHONY: force
NP := $(shell procsel $(PETSCFEM_PROCTABLE))
run: $(CASE_NAME).depl
	-$(MPIRUN) -np 1 -machinefile machi.dat $(PROG) \
			-case $<  > output.$(case_name).tmp
	if [ $$oct_check = yes ] ; then \
		 $(OCTAVE) check.m > verif.$(case_name).tmp ; fi

genload.depl: force

constraint:
	$(MAKE) case_name=$@ run

scalar:
	$(MAKE) case_name=$@ run

diag:
	$(MAKE) case_name=$@ run

full:
	$(MAKE) case_name=$@ run

fstack_global:
	$(MAKE) case_name=$@ oct_check=no run

local_clean::
	-rm -rf mat.output 

#s
