#__INSERT_LICENSE__
# $Id: Makefile,v 1.30 2003/03/10 03:00:21 mstorti Exp $ 

.PHONY: force all run tests data viscous inviscid sphere

default: nsid

tests: 

nsid:
	$(MAKE) CASE_NAME=nsid run

sphere:
	$(MAKE) CASE_NAME=sphere run

inviscid:
	$(MAKE) CASE_NAME=ext NP=1 run > inviscid.output.tmp 2>&1 

plate:
	$(MAKE) CASE_NAME=plate run

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
# Currently the coupled version doesn't work in parallel
# ifeq ($(NP),)
# NP := $(shell procsel ./proctable)
# endif
NP := 1
run: $(CASE_NAME).depl
	-$(MPIRUN) -np $(NP) -machinefile machi.dat $(PROG) -case $< 

$(CASE_NAME).depl: force

local_clean::
	-rm -rf mat.output *.fifo

dx:
	eperl -P viewer.dx > viewer.dx.tmp
#s
