# $Id: Makefile,v 1.42 2003/02/11 21:07:43 mstorti Exp $ 

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
PROG := $(NS_O)

CASE := $(CASE_NAME)
################################################################
ifeq (1,0)
NP = $(shell procsel -n $(PETSCFEM_PROCTABLE))
# NP = $(shell procsel -n ./proctable.geronimo)
NL = -nolocal
MCH = -machinefile ./machi.dat
endif

ifeq (1,1)
NP = 1
MCH =
endif

export NP

PROG := $(NS)

#w Runs program
$(CASE).depl: force
run: $(CASE).depl
#	$(MAKE) mesh
	-${MPIRUN} $(NL) -np $(NP) $(MCH) $(PROG) -case $(CASE).depl
#	-${MPIRUN} -np $(NP) $(MCH) $(PROG) -case $(CASE).depl

mesh: cubcav.depl
	$(OCTAVE) mkmesh.m

local_clean::
	-chmod +w *.depl
	-rm *.depl *.sal tmp_file*.tmp &> /dev/null

hexasplit.bin: hexasplit.o
	${CXX_CLINKER} -o $@ $<  $(LDFLAGS_NO_PROG_LIB)

nettest = $(PETSCFEM_DIR)/tools/nettest.bin

.PHONY: nettest
nettest:
	NP=`procsel -n ./proctable.geronimo` ;  \
	$(MPIRUN) -np $$NP -machinefile machi.dat $(nettest)

dx:
	eperl -P viewer.dx > viewer.dx.tmp
