# $Id: Makefile,v 1.20 2002/07/28 00:21:56 mstorti Exp $ 

.PHONY: all run  compile lclean laplace adv nso laplace_tests force \
        ctest gprof

all: run

#default case
# case := laplace
case := iisd_sbp1_uj
export case

tests: test_re1000 test_iisd

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

NS = $(PETSCFEM_DIR)/applications/ns/ns_O.bin

LOC_PROG := ns
#CASE_NAME := qharm
CASE_NAME := cubcav
PROG := $(NS)

CASE := $(CASE_NAME)
################################################################
NP = 1

PROG := $(NS)
MCH = -machinefile ./machi.dat

#w Runs program
$(CASE).depl: force
run: $(CASE).depl
#	$(MAKE) mesh
#	-${MPIRUN} -nolocal -np $(NP) $(MCH) $(PROG) -case $(CASE).depl
	-${MPIRUN} -np $(NP) $(MCH) $(PROG) -case $(CASE).depl

mesh: cubcav.depl
	octave -q mkmesh.m

local_clean::
	-chmod +w *.depl
	-rm *.depl *.sal tmp_file*.tmp &> /dev/null

