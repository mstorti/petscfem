#__INSERT_LICENSE__
# $Id: Makefile,v 1.3 2001/05/25 21:29:54 mstorti Exp $ 

.PHONY: force all run tests data

all: run

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include ../Makefile.base.test
SRCS = 

#s Test/Advdif tests
tests: test_case_1 test_case_2a test_case_2b

.PHONY: data force
NP := $(shell procsel $(PETSCFEM_PROCTABLE))
test_case_%: 
	$(MAKE) EPERLFLAGS='-d case=$*' data
	-$(MPIRUN) -np 1 -machinefile machi.dat $(ADVDIF) \
			-case genload.depl > output.tmp
	mv save.state.tmp save.case_$*.tmp

data: genload.depl

genload.depl: force

local_clean::
	-rm -rf mat.output 

auxi:
	if [ ! -L advdif ] ; \
	 then ln -sf $(PETSCFEM_DIR)/applications/advdif/advdif.bin advdif ; \
	fi
	echo file advdif >.gdbinit
	echo set args -case genload.depl >>.gdbinit
#s
