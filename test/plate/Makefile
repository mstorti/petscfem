#__INSERT_LICENSE__
# $Id: Makefile,v 1.7 2002/12/30 03:06:31 mstorti Exp $ 

.PHONY: force all run tests data viscous inviscid

all: run

tests: 

inviscid:
	$(MAKE) CASE_NAME=ext run

viscous:
	if [ -e visc2inv.fifo ] ; then rm visc2inv.fifo ; fi
	mkfifo visc2inv.fifo
	if [ -e inv2visc.fifo ] ; then rm inv2visc.fifo ; fi
	mkfifo inv2visc.fifo
	$(MAKE) CASE_NAME=cylin run

SRCS = inviscid.cpp
PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

LOC_PROG := ns
PROG := $(NS_O)

.PHONY: force
NP := $(shell procsel $(PETSCFEM_PROCTABLE))
run: $(CASE_NAME).depl
	-$(MPIRUN) -np 1 -machinefile machi.dat $(PROG) -case $< 

$(CASE_NAME).depl: force

local_clean::
	-rm -rf mat.output visc2inv inv2visc

dx:
	eperl -P viewer.dx > viewero.dx
#s
