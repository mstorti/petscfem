#__INSERT_LICENSE__
# $Id: Makefile,v 1.3 2002/12/22 06:20:44 mstorti Exp $ 

.PHONY: force all run tests data scalar diag full

default: scalar

tests: scalar diag full

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test
SRCS = 

LOC_PROG := ns
CASE_NAME := genload
PROG := $(NS_G)

.PHONY: force
NP := $(shell procsel $(PETSCFEM_PROCTABLE))
run: $(CASE_NAME).depl
	-$(MPIRUN) -np 1 -machinefile machi.dat $(PROG) \
			-case $< #  > output.$(CASE_NAME).tmp
#	octave -qH check.m  # > verif.$(case_name).tmp

genload.depl: force

scalar:
	$(MAKE) case_name=$@ run

diag:
	$(MAKE) case_name=$@ run

full:
	$(MAKE) case_name=$@ run


local_clean::
	-rm -rf mat.output 

#s
