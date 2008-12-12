#__INSERT_LICENSE__

.PHONY: force all run tests data scalar diag full constraint default \
	test_fstack fstack_global fstack_nodes fstack_nodes2 qconst

default: qconst

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
# MCH := -machinefile machi.dat
run: $(CASE_NAME).depl
	-$(MPIEXEC) $(MCH) -n $(NP) $(PROG) \
			-case $< $(call redir,output.$(case_name).tmp)
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

qconst:
	$(MAKE) case_name=$@ oct_check=no run

fstack_nodes fstack_global fstack_nodes2:
	$(MAKE) case_name=$@ oct_check=no NP=2 run

local_clean::
	-rm -rf mat.output 

#s
