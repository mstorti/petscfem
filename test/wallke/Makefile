#__INSERT_LICENSE__
# $Id: Makefile,v 1.3 2001/06/04 22:41:12 mstorti Exp $ 

.PHONY: force all run tests data

all: run

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include ../Makefile.base.test
SRCS = 

ifeq ($(case),"test_case_a")
NAME = nonlr
else
NAME = wallke
endif

PROG = $(NS)
LOC_PROG = ns

tests: test_case_a

.PHONY: data force
NP := $(shell procsel $(PETSCFEM_PROCTABLE))
test_case_%_data: 
	-if [ ! "$(cont)" = "y" ] ; then rm -f *.tmp ; fi 
	$(MAKE) EPERLFLAGS='-d case=$* -d cont=$(cont)' data
	octave -q mkmesh.m

data: $(NAME).depl force
$(NAME).depl: force

#data.m.tmp: $(NAME).epl
#	$(MAKE) $(NAME).depl

test_case_%: test_case_%_data
	-$(MPIRUN) -np 1 -machinefile machi.dat $(PROG) \
			-case $(NAME).depl #> output.case_$*.tmp
	cp save.state save.state.tmp

local_clean::
	-rm -rf mat.output system.dat save.state data.m

auxi:
	if [ ! -L advdif ] ; \
	 then ln -sf $(PROG) $(LOC_PROG) ; \
	fi
	echo file ns >.gdbinit
	echo set args -case $(NAME).depl >>.gdbinit
#s
