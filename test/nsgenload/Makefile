#__INSERT_LICENSE__
# $Id: Makefile,v 1.10 2003/02/14 00:39:02 mstorti Exp $ 

.PHONY: force all run tests data scalar diag full constraint default \
	test_fstack fstack_global fstack_nodes fstack_nodes2

default: constraint

tests: scalar diag full constraint test_fstack

test_fstack: fstack_global fstack_nodes fstack_nodes2 

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test
SRCS = 

LOC_PROG := ns
CASE_NAME := genload
PROG := $(NS)

oct_check := yes
export oct_check

ifeq ($(NP),)
NP := 1
endif

NPP := $(shell procsel $(PETSCFEM_PROCTABLE))
run: $(CASE_NAME).depl
	-$(MPIRUN) -np $(NP) -machinefile machi.dat $(PROG) \
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

fstack_nodes fstack_global fstack_nodes2:
	$(MAKE) case_name=$@ oct_check=no NP=2 run

local_clean::
	-rm -rf mat.output 

#s
