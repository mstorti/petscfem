#__INSERT_LICENSE__
# $Id: Makefile,v 1.6 2002/12/29 22:27:24 mstorti Exp $ 

.PHONY: force all run tests data scalar diag full cylin

all: run

tests: 

SRCS = inviscid.cpp
PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

LOC_PROG := ns
# CASE_NAME := plate
CASE_NAME := cylin
PROG := $(NS_O)

.PHONY: force
NP := $(shell procsel $(PETSCFEM_PROCTABLE))
run: $(CASE_NAME).depl
	if [ -e visc2inv ] ; then rm visc2inv ; fi
	mkfifo visc2inv
	if [ -e inv2visc ] ; then rm inv2visc ; fi
	mkfifo inv2visc
	-$(MPIRUN) -np 1 -machinefile machi.dat $(PROG) -case $< 

$(CASE_NAME).depl: force

local_clean::
	-rm -rf mat.output visc2inv inv2visc

dx:
	eperl -P viewer.dx > viewero.dx
#s
