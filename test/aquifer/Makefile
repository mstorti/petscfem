# $Id: Makefile,v 1.1 2002/01/15 21:35:40 mstorti Exp $ 

.PHONY:  force

tests: 

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

LOC_PROG := advdif
CASE_NAME := aqui
PROG := $(ADVDIF)

CASE = $(CASE_NAME)
################################################################
NP = 1

PROG := $(ADVDIF)
MCH = -machinefile ./machi.dat

NP := $(shell procsel $(PETSCFEM_PROCTABLE))
ifeq ($(np),)
np := 1
endif

#w Runs program
run:  data
	-${MPIRUN} -np $(np) $(MCH) $(PROG) -case $(CASE).depl 

data: force
	$(MAKE) $(CASE).depl
	octave -q mkmesh.m

# EPERLFLAGS = -d case=$$case
$(CASE).depl: force
#	echo in target \$(CASE).depl : CASE $(CASE)

local_clean::
	-chmod +w *.depl
	-rm *.depl *.sal tmp_file*.tmp &> /dev/null
