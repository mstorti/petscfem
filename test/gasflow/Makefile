#__INSERT_LICENSE__
# $Id: Makefile,v 1.5 2003/01/16 02:08:22 mstorti Exp $ 

.PHONY: force all run tests

all: run

tests: 

gasflow:
	$(MAKE) CASE_NAME=gasflow_test_2d run

SRCS = 
PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

ifeq ($(CASE_NAME),) 
CASE_NAME := nozzle
# CASE_NAME := gasflow_test_2d
endif
LOC_PROG := advdif
PROG := ../../applications/advdif/advdif_O.bin

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
