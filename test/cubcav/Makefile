# $Id: Makefile,v 1.49 2003/08/25 22:20:21 mstorti Exp $ 

.PHONY: all run  compile lclean laplace adv nso laplace_tests force \
        ctest gprof

default: 
	$(MAKE) CASE=cubcav run

test_isp: 

#default case
# case := laplace
# case := iisd_sbp1_uj
export case

tests: test_re1000 test_iisd

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

LOC_PROG := ns
#CASE_NAME := qharm
CASE_NAME := cubcav
PROG := $(NS)

CASE := $(CASE_NAME)
################################################################
ifeq (1,1)
# NP = $(shell procsel -n $(PETSCFEM_PROCTABLE))
NP = $(shell procsel -n ./proctable)
NL = -nolocal
MCH = -machinefile ./machi.dat
endif

ifeq (1,0)
NP = 1
MCH =
endif

NP:=1

export NP

PROG := $(NS)

#w Runs program
run: 
	make-depl cubcav.epl cubcav.depl -d case=ns
	-${MPIRUN} $(NL) -np $(NP) $(MCH) $(PROG) -case $(CASE).depl

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

dxo:
	$(EPERL) -P viewer.dx > viewer.dx.tmp

dx:
	$(EPERL) -P -d dx=1 cubcav.epl > cubcav_dx.depl
	${MPIRUN} -np 1 $(PROG) -case cubcav_dx.depl

k:
	-killall ns_O.bin job.pl
