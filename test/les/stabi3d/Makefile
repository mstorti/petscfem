#__INSERT_LICENSE__
# $Id: Makefile,v 1.3 2002/11/10 14:26:43 mstorti Exp $ 

.PHONY: force tangent stabi data auxi

run: stabi

CASE := stabi
PETSCFEM_DIR = ../../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

stabi: 	
	$(MAKE) $(CASE).depl
	octave -q < mkmesh.m 
	-$(MPIRUN) -np 1 $(PROG) -case stabi.depl # >stabi.out.tmp

LOC_PROG := ns
CASE_NAME := stabi
PROG := $(PETSCFEM_DIR)/applications/ns/ns_O.bin
