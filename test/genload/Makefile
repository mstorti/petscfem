#__INSERT_LICENSE__

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
test_case_%: 
	$(MAKE) EPERLFLAGS='-d case=$*' data
	-$(MPIEXEC) -n 1 $(ADVDIF) \
			-case genload.depl </dev/null > output.case_$*.tmp
	if [ -e save.state.tmp ] ; then mv save.state.tmp save.case_$*.tmp ; fi

data: genload.depl

genload.depl: force

local_clean::
	-rm -rf mat.output 

LOC_PROG := advdif
#CASE_NAME := qharm
CASE_NAME := genload
PROG := $(ADVDIF)

link:
	ln -sf $(PROG) .

auxi_o:
	if [ ! -L advdif ] ; \
	 then ln -sf $(PETSCFEM_DIR)/applications/advdif/advdif_g.bin advdif ; \
	fi
	echo file advdif_g >.gdbinit
	echo set args -case genload.depl >>.gdbinit
#s
