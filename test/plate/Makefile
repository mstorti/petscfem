#__INSERT_LICENSE__
# $Id: Makefile,v 1.20 2003/02/24 23:19:34 mstorti Exp $ 

.PHONY: force all run tests data viscous inviscid sphere

default: wave

tests: 

sphere:
	$(MAKE) CASE_NAME=sphere run

inviscid:
	$(MAKE) CASE_NAME=ext run > inviscid.output.tmp 2>&1 

idx:
	$(MAKE) CASE_NAME=ext run

viscous: inviscid.efn
	if [ -e visc2inv.fifo ] ; then rm visc2inv.fifo ; fi
	mkfifo visc2inv.fifo
	if [ -e inv2visc.fifo ] ; then rm inv2visc.fifo ; fi
	mkfifo inv2visc.fifo
	$(MAKE) CASE_NAME=cylin run

wave:
	$(MAKE) CASE_NAME=wave run

inviscid.efn: EFN_LDFLAGS = ./fifo.o
inviscid.efn: fifo.o
fifo.o: fifo.cpp

SRCS = inviscid.cpp
PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

LOC_PROG := ns
PROG := $(NS)

.PHONY: force
NP := $(shell procsel $(PETSCFEM_PROCTABLE))
run: $(CASE_NAME).depl
	-$(MPIRUN) -np 1 -machinefile machi.dat $(PROG) -case $< 

$(CASE_NAME).depl: force

local_clean::
	-rm -rf mat.output *.fifo

dx:
	eperl -P viewer.dx > viewer.dx.tmp
#s
