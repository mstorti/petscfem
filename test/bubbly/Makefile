# $Id: Makefile,v 1.5 2002/02/27 15:05:13 mstorti Exp $ 

.PHONY: all run  compile lclean laplace adv nso laplace_tests force \
        ctest gprof

tests: test_sqcav

all: run

# PETSCFEM_DIR = $(HOME)/PETSC/petscfem
PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

LOC_PROG := advdif_g
# CASE_NAME := bubbly
CASE_NAME := couette
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
	$(MAKE) $(CASE).depl
	octave -q mkcouette.m

# EPERLFLAGS = -d case=$$case
$(CASE).depl: force
#	echo in target \$(CASE).depl : CASE $(CASE)

local_clean::
	-chmod +w *.depl
	-rm *.depl *.sal tmp_file*.tmp &> /dev/null
