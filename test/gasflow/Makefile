#__INSERT_LICENSE__
# $Id: Makefile,v 1.10 2003/01/21 16:09:16 mstorti Exp $ 

.PHONY: force all run tests default strip gasflow nozzle  vtubei

default: vtubei

tests: 

gasflow:
	$(MAKE) CASE_NAME=gasflow_test_2d run

nozzle:
	$(MAKE) CASE_NAME=$@ run

vtube:
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
PROG := $(ADVDIF)
endif

# NP := $(shell procsel $(PETSCFEM_PROCTABLE))
run: $(CASE_NAME).depl
	-$(MPIRUN) -np 1 -machinefile machi.dat $(PROG) -case $< |  \
		 tee $(CASE_NAME).output.tmp

$(CASE_NAME).depl: force

local_clean::
	-rm -rf mat.output visc2inv inv2visc

dx:
	eperl -P viewer.dx > viewero.dx
#s
