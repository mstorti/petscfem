#__INSERT_LICENSE__
# $Id: Makefile,v 1.2 2001/05/31 02:32:27 mstorti Exp $ 

.PHONY: force all run tests data

all: run

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include ../Makefile.base.test
SRCS = 

NAME = nonlr
PROG = $(NS)
LOC_PROG = ns

tests: test_case_a

.PHONY: data force
NP := $(shell procsel $(PETSCFEM_PROCTABLE))
test_case_%: 
	$(MAKE) EPERLFLAGS='-d case=$*' data
	-$(MPIRUN) -np 1 -machinefile machi.dat $(PROG) \
			-case $(NAME).depl #> output.case_$*.tmp
#	if [ -e save.state.tmp ] ; then mv save.state.tmp save.case_$*.tmp ; fi

data: $(NAME).depl

$(NAME).depl: force

local_clean::
	-rm -rf mat.output system.dat 

auxi:
	if [ ! -L advdif ] ; \
	 then ln -sf $(PROG) $(LOC_PROG) ; \
	fi
	echo file ns >.gdbinit
	echo set args -case $(NAME).depl >>.gdbinit
#s
