#__INSERT_LICENSE__
# $Id: Makefile,v 1.7 2003/03/20 10:14:08 mstorti Exp $ 

.PHONY: force all run tests data viscous inviscid sphere

default: spiller

tests: 

spiller:
	$(MAKE) CASE_NAME=spiller run

spillway_mesh:
	initia=0 octave -qH spiller.m

SRCS = inviscid.cpp
PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

LOC_PROG := ns
PROG := $(NS)

.PHONY: force
# NP := $(shell procsel ./proctable.tmp)
NP := 1
run: $(CASE_NAME).depl alehook.efn
	-$(MPIRUN) -np $(NP) -machinefile machi.dat $(PROG) -case $< 

$(CASE_NAME).depl: force

local_clean::
	-rm -rf mat.output *.fifo

#s
