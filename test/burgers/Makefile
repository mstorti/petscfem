#__INSERT_LICENSE__
# $Id: Makefile,v 1.2 2002/07/25 22:35:10 mstorti Exp $ 

.PHONY: force all run tests data

all: run

test_%:
	$(MAKE) CASE=$* run

tests: test__d05 test__d01 test__d001

test__d01: test__d01_wf test__d01_nwf 
	CASE=_d01 $(OCTAVE) -q verif.m > $@.out.tmp

test__d05: test__d05_wf test__d05_nwf 
	CASE=_d05 $(OCTAVE) -q verif.m > $@.out.tmp

test__d001: test__d001_wf test__d001_nwf 
	CASE=_d001 $(OCTAVE) -q verif.m > $@.out.tmp

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include ../Makefile.base.test
SRCS = 

#s Test/Advdif tests
# ADVDIF    = $(PETSCFEM_DIR)/applications/advdif/advdif.bin 

OF = burgers.nod.tmp burgers.con.tmp burgers.ini.tmp
$(OF): mmesh.m
	$(OCTAVE) -q mmesh.m

data: burgers.depl $(OF) 
burgers.depl: burgers.epl data.m force

################################################################
# TESTS FOR THE NEW FLUX FUNCTIONS...

NP := $(shell procsel $(PETSCFEM_PROCTABLE))
run: data
	-$(MPIRUN) -np 1 -machinefile machi.dat $(ADVDIF) \
			-case burgers.depl > burgers.out
	mv -f save.state.tmp save.state.$(CASE).tmp

local_clean::

#s
