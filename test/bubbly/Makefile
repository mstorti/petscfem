# $Id: Makefile,v 1.2 2002/02/17 03:59:07 mstorti Exp $ 

.PHONY: all run  compile lclean laplace adv nso laplace_tests force \
        ctest gprof

tests: test_sqcav

all: run

# PETSCFEM_DIR = $(HOME)/PETSC/petscfem
PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

LOC_PROG := advdif
#CASE_NAME := qharm
CASE_NAME := bubbly
PROG := $(ADVDIF)

CASE := $(CASE_NAME)

################################################################
NP = 1

MCH = -machinefile ./machi.dat

NP := $(shell procsel $(PETSCFEM_PROCTABLE))
ifeq ($(np),)
np := 1
endif

#w Runs program
run:  data
	-${MPIRUN} -np $(np) $(MCH) $(PROG) \
			-case $(CASE).depl # > output.out

data: force
	echo in target data: CASE = $(CASE)
	$(MAKE) $(CASE).depl
	echo "in target data : CASE $(CASE)"
	octave -q mkmesh.m

# EPERLFLAGS = -d case=$$case
$(CASE).depl: force
#	echo in target \$(CASE).depl : CASE $(CASE)

local_clean::
	-chmod +w *.depl
	-rm *.depl *.sal tmp_file*.tmp &> /dev/null
