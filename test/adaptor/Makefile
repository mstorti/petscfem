#__INSERT_LICENSE__
# $Id: Makefile,v 1.1.2.1 2001/10/29 14:34:41 mstorti Exp $ 

.PHONY: force all run tests data

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base

NSO := $(PETSCFEM_DIR)/../petscfem-beta-1.98/applications/ns/ns.bin

NPP := $(shell procsel $(PETSCFEM_PROCTABLE))

PROG := $(NS)
run: 
	octave < mkmesh.m
#	-$(MPIRUN) -np 1 -machinefile machi.dat $(PROG) \
#			-case elast.dat

local_clean::
	-rm -rf mat.output system.dat save.state data.m

LOC_PROG := ns
CASE_NAME := elast
ifeq ($(PROG),)
PROG := $(NS)
endif

#s
