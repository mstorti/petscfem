#__INSERT_LICENSE__
# $Id: Makefile,v 1.9 2003/01/08 15:54:26 mstorti Exp $ 

.PHONY: force all run tests data

all: test_case_poiseuille

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test
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
	$(OCTAVE) mkmesh.m

data: $(NAME).depl force
$(NAME).depl: force

#data.m.tmp: $(NAME).epl
#	$(MAKE) $(NAME).depl

test_case_%: test_case_%_data
	if [ -e save.state ] ; then cp save.state save.state.tmp ; fi
	-$(MPIEXEC) -n 1 -machinefile machi.dat $(PROG) \
			-case $(NAME).depl | tee output.case_$*.tmp

local_clean::
	-rm -rf mat.output system.dat save.state data.m

auxi_o:
	if [ ! -L advdif ] ; \
	 then ln -sf $(PROG) $(LOC_PROG) ; \
	fi
	echo file ns >.gdbinit
	echo set args -case $(NAME).depl >>.gdbinit
#s
