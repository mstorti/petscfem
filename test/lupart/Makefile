#__INSERT_LICENSE__
# $Id: Makefile,v 1.4 2001/07/18 22:44:06 mstorti Exp $ 

.PHONY: force all run tests data

all: test_case_poiseuille
dat: test_case_poiseuille_data

#PETSCFEM_DIR = ../..
PETSCFEM_DIR = $(HOME)/PETSC/petscfem-beta-1.82
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test
SRCS = 
NAME = wallke

PROG = $(NS)
LOC_PROG = ns

tests: test_case_a

.PHONY: data force
test_case_%_data: 
	-if [ ! "$(cont)" = "y" ] ; then rm -f *.tmp ; fi 
	$(MAKE) EPERLFLAGS='-d case=$* -d cont=$(cont)' data
	octave -q mkmesh.m

data: $(NAME).depl force
$(NAME).depl: force

NPP := $(shell procsel $(PETSCFEM_PROCTABLE))
NPP := 1
ifeq ($(NP),)
NP := $(NPP)
endif
test_case_%: test_case_%_data
	if [ -e save.state ] ; then cp save.state save.state.tmp ; fi
	-$(MPIRUN) -np $(NP) -machinefile machi.dat $(PROG) \
			-case $(NAME).depl | tee output.case_$*.tmp

local_clean::
	-rm -rf mat.output system.dat save.state data.m

auxi:
	if [ ! -L advdif ] ; \
	 then ln -sf $(PROG) $(LOC_PROG) ; \
	fi
	echo file ns >.gdbinit
	echo set args -case $(NAME).depl >>.gdbinit
#s
