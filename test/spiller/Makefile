#__INSERT_LICENSE__
# $Id: Makefile,v 1.3 2003/03/17 22:49:31 mstorti Exp $ 

.PHONY: force all run tests data viscous inviscid sphere

default: spiller

tests: 

spiller:
	$(MAKE) CASE_NAME=spiller run

SRCS = inviscid.cpp
PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

LOC_PROG := ns
PROG := $(NS)

.PHONY: force
# NP := $(shell procsel ./proctable.tmp)
NP := 1
run: $(CASE_NAME).depl
	-$(MPIRUN) -np $(NP) -machinefile machi.dat $(PROG) -case $< 

$(CASE_NAME).depl: force

local_clean::
	-rm -rf mat.output *.fifo

#s
