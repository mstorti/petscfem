# $Id: Makefile,v 1.1 2005/01/07 23:39:12 mstorti Exp $ 

.PHONY: force all run tests default strip gasflow nozzle  vtubei

default: run

SRCS = 
PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

ifeq ($(CASE_NAME),) 
CASE_NAME := gfplate
endif
LOC_PROG := advdif
ifeq ($(PROG),)
PROG := $(ADVDIF)
NL :=
NP := 1
endif

run: 
	make-depl $(CASE_NAME).epl $(CASE_NAME).depl
	-$(MPIRUN) $(NL) -np $(NP) -machinefile		\
	machi.dat $(PROG) -case $(CASE_NAME).depl # > $(CASE_NAME).log

link:
	$(MAKE) -C ../../ $(LOC_PROG)
	ln -f $(PROG) $(LOC_PROG)

#s
