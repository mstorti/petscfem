# $Id: Makefile,v 1.2 2003/11/10 21:29:50 mstorti Exp $ 

.PHONY: run default

default: run

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

LOC_PROG := advdif
CASE := pool
PROG := $(ADVDIF)

###--------<*>--------<*>--------<*>--------<*>--------<*>--------<*>--------<*> 
NP = 1

MCH = -machinefile ./machi.dat

NP := $(shell procsel $(PETSCFEM_PROCTABLE))
ifeq ($(np),)
np := 1
endif

#w Runs program
run:
	make-depl $(CASE).epl $(CASE).depl
	octave -qH mkmesh.m
	-${MPIRUN} -np $(np) $(MCH) $(PROG) -case $(CASE).depl
