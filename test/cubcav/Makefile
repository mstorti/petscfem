# $Id: Makefile,v 1.29 2002/09/12 21:41:05 mstorti Exp $ 

.PHONY: all run  compile lclean laplace adv nso laplace_tests force \
        ctest gprof

all: run

#default case
# case := laplace
case := iisd_sbp1_uj
export case

tests: test_re1000 test_iisd

PETSCFEM_DIR = ../..
# PETSCFEM_DIR = $(HOME)/PETSC/petscfem-beta-3.03.pl5
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

NS = $(PETSCFEM_DIR)/applications/ns/ns_O.bin

LOC_PROG := ns
#CASE_NAME := qharm
CASE_NAME := cubcav
PROG := $(NS)

CASE := $(CASE_NAME)
################################################################
ifeq (1,1)
NP = $(shell procsel -n ./proctable.geronimo)
NL = -nolocal
endif

ifeq (1,0)
NP = $(shell procsel ./proctable)
NP = 2
endif

# NP = 2
export NP

PROG := $(NS)
MCH = -machinefile ./machi.dat

#w Runs program
$(CASE).depl: force
run: $(CASE).depl
#	$(MAKE) mesh
	-${MPIRUN} $(NL) -np $(NP) $(MCH) $(PROG) -case $(CASE).depl
#	-${MPIRUN} -np $(NP) $(MCH) $(PROG) -case $(CASE).depl

mesh: cubcav.depl
	octave -q mkmesh.m

local_clean::
	-chmod +w *.depl
	-rm *.depl *.sal tmp_file*.tmp &> /dev/null

hexasplit.bin: hexasplit.o
	${CXX_CLINKER} -o $@ $<  $(LDFLAGS_NO_PROG_LIB)

