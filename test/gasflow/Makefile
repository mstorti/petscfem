# $Id: Makefile,v 1.18 2003/01/27 00:37:11 mstorti Exp $ 

.PHONY: force all run tests default strip gasflow nozzle  vtubei

default: vtube

tests: 

gasflow:
	$(MAKE) CASE_NAME=gasflow_test_2d run

nozzle:
	$(MAKE) CASE_NAME=$@ run

vtube: gasflow.efn
	$(MAKE) CASE_NAME=$@ run

vtubei:
	$(MAKE) CASE_NAME=vtube incomp=1 PROG=$(NS_O) run

strip:
	$(MAKE) CASE_NAME=$@ run

SRCS = 
PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

ifeq ($(CASE_NAME),) 
CASE_NAME := strip
endif
LOC_PROG := advdif
ifeq ($(PROG),)
PROG := $(ADVDIF_O)
endif

# NP := $(shell procsel -n $(PETSCFEM_PROCTABLE))
include ../les/stabi3d/Makefile.parproc
run: $(CASE_NAME).depl
	-$(MPIRUN) -nolocal -np $(NP) -machinefile machi.dat $(PROG) -case $< 

$(CASE_NAME).depl: force

local_clean::
	-rm -rf mat.output visc2inv inv2visc

dx:
	eperl -P viewer.dx > viewer_tmp.dx

gasflow.efn: gasflow.cpp
#s
