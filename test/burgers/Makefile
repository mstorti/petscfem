#__INSERT_LICENSE__
# $Id: Makefile,v 1.8 2003/11/15 16:04:10 mstorti Exp $ 

.PHONY: force all run tests data

default: test__d01_wf

all: run

test_%:
	$(MAKE) CASE=$* run

#w tests: test__d05 test__d01 test__d001
tests: test__d01 test__d001

#w Burger test with viscosity 0.1
test__d01: test__d01_wf test__d01_nwf 
	CASE=_d01 $(OCTAVE) verif.m > $@.out.tmp

#w Burger test with viscosity 0.5
test__d05: test__d05_wf test__d05_nwf 
	CASE=_d05 $(OCTAVE) verif.m > $@.out.tmp

#w Burger test with viscosity 0.01
test__d001: test__d001_wf test__d001_nwf 
	CASE=_d001 $(OCTAVE) verif.m > $@.out.tmp

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include ../Makefile.base.test
SRCS = 

#s Test/Advdif tests
# ADVDIF    = $(PETSCFEM_DIR)/applications/advdif/advdif.bin 

OF = burgers.nod.tmp burgers.con.tmp burgers.ini.tmp
$(OF): mmesh.m data.m
	$(OCTAVE) mmesh.m

data: burgers.depl $(OF) 
burgers.depl: burgers.epl data.m force

################################################################
# TESTS FOR THE NEW FLUX FUNCTIONS...

NP := $(shell procsel $(PETSCFEM_PROCTABLE))
run: data
	-$(MPIRUN) -np 1 -machinefile machi.dat $(ADVDIF) \
			-case burgers.depl > burgers.$(CASE).out.tmp
	mv -f save.state.tmp save.state.$(CASE).tmp

local_clean::

PROG := $(ADVDIF)
LOC_PROG := advdif

#w Makes local link for running with debugger
link:
	$(MAKE) -C $(PETSCFEM_DIR) $(LOC_PROG)
	ln -f $(PROG) $(LOC_PROG)

#s
