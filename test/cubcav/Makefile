# $Id: Makefile,v 1.2.2.1 2002/01/14 02:17:14 mstorti Exp $ 

.PHONY: all run  compile lclean laplace adv nso laplace_tests force \
        ctest gprof

all: run

tests: test_re1000 test_iisd

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

NS = $(PETSCFEM_DIR)/applications/ns/ns.bin

LOC_PROG := ns
#CASE_NAME := qharm
CASE_NAME := cubcav
PROG := $(NS)

CASE := $(CASE_NAME)
################################################################
NP = 1

PROG := $(NS)
MCH = -machinefile ./machi.dat

NP := $(shell procsel $(PETSCFEM_PROCTABLE))

#w Runs program
case := iisd_sbp1_uj
export case
$(CASE).depl: force
run: $(CASE).depl
	octave -q mkmesh.m
	-${MPIRUN} -np 1 $(MCH) $(PROG) -case $(CASE).depl

local_clean::
	-chmod +w *.depl
	-rm *.depl *.sal tmp_file*.tmp &> /dev/null
