#__INSERT_LICENSE__
# $Id: Makefile,v 1.2 2003/03/21 20:27:50 mstorti Exp $ 

.PHONY: force all run tests data viscous inviscid sphere

default: spillway

tests: 

spillway:
	$(MAKE) CASE_NAME=spillway run

spillway_mesh:
	initia=0 octave -qH spillway.m

SRCS = inviscid.cpp
PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

LOC_PROG := ns
PROG := $(NS)

.PHONY: force
# NP := $(shell procsel $(PETSCFEM_PROCTABLE))
# NP := 1
NP := $(shell procsel ./proctable.tmp)
run: $(CASE_NAME).depl alehook.efn
	-$(MPIRUN) -np $(NP) -machinefile machi.dat $(PROG) -case $< 

$(CASE_NAME).depl: force

local_clean::
	-rm -rf mat.output *.fifo

#s
