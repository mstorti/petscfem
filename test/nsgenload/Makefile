#__INSERT_LICENSE__
# $Id: Makefile,v 1.8 2003/01/08 15:54:26 mstorti Exp $ 

.PHONY: force all run tests data scalar diag full constraint default

default: constraint

tests: scalar diag full constraint

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test
SRCS = 

LOC_PROG := ns
CASE_NAME := genload
PROG := $(NS)

.PHONY: force
NP := $(shell procsel $(PETSCFEM_PROCTABLE))
run: $(CASE_NAME).depl
	-$(MPIRUN) -np 1 -machinefile machi.dat $(PROG) \
			-case $<  > output.$(CASE_NAME).tmp
	$(OCTAVE) check.m > verif.$(case_name).tmp

genload.depl: force

constraint:
	$(MAKE) case_name=$@ run

scalar:
	$(MAKE) case_name=$@ run

diag:
	$(MAKE) case_name=$@ run

full:
	$(MAKE) case_name=$@ run


local_clean::
	-rm -rf mat.output 

#s
