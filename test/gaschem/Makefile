# $Id: Makefile,v 1.3 2003/11/11 02:13:00 mstorti Exp $ 

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

#w Makes local link for running with debugger
link:
	$(MAKE) -C $(PETSCFEM_DIR) $(LOC_PROG)
	ln -f $(PROG) $(LOC_PROG)
