#__INSERT_LICENSE__
# $Id: Makefile,v 1.4 2002/11/11 12:50:19 mstorti Exp $ 

.PHONY: force tangent stabi data auxi

run: stabi

CASE := stabi
PETSCFEM_DIR = ../../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

NP := $(shell procsel -n $(PETSCFEM_PROCTABLE))

stabi: 	
	$(MAKE) $(CASE).depl
	octave -q < mkmesh.m 
	-$(MPIRUN) -machinefile ./machi.dat \
		-np $(NP) -nolocal $(PROG) -case stabi.depl # >stabi.out.tmp

LOC_PROG := ns
CASE_NAME := stabi
PROG := $(PETSCFEM_DIR)/applications/ns/ns_O.bin
