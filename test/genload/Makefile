#__INSERT_LICENSE__
# $Id: Makefile,v 1.6 2001/07/28 20:04:07 mstorti Exp $ 

.PHONY: force all run tests data

all: run

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include ../Makefile.base.test
SRCS = 

#s Test/Advdif tests
tests: test_case_1 test_case_2a test_case_2b test_case_fstack \
	test_case_genl1d test_case_genl1dh

.PHONY: data force
NP := $(shell procsel $(PETSCFEM_PROCTABLE))
test_case_%: 
	$(MAKE) EPERLFLAGS='-d case=$*' data
	-$(MPIRUN) -np 1 -machinefile machi.dat $(ADVDIF) \
			-case genload.depl > output.case_$*.tmp
	if [ -e save.state.tmp ] ; then mv save.state.tmp save.case_$*.tmp ; fi

data: genload.depl

genload.depl: force

local_clean::
	-rm -rf mat.output 

auxi_o:
	if [ ! -L advdif ] ; \
	 then ln -sf $(PETSCFEM_DIR)/applications/advdif/advdif.bin advdif ; \
	fi
	echo file advdif >.gdbinit
	echo set args -case genload.depl >>.gdbinit
#s
