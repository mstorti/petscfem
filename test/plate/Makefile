#__INSERT_LICENSE__
# $Id: Makefile,v 1.2 2002/12/19 18:51:25 mstorti Exp $ 

.PHONY: force all run tests data scalar diag full

all: run

tests: 

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test
SRCS = 

LOC_PROG := ns
CASE_NAME := plate
PROG := $(NS_G)

.PHONY: force
NP := $(shell procsel $(PETSCFEM_PROCTABLE))
run: $(CASE_NAME).depl
	-$(MPIRUN) -np 1 -machinefile machi.dat $(PROG) \
			-case $< #  > output.$(CASE_NAME).tmp

$(CASE_NAME).depl: force

local_clean::
	-rm -rf mat.output 

dx:
	eperl -P viewer.dx > viewero.dx
#s
