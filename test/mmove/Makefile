#__INSERT_LICENSE__
# $Id: Makefile,v 1.1.2.1 2001/12/20 02:32:01 mstorti Exp $ 

.PHONY: force all run tests data

PETSCFEM_DIR = ../..
#PETSCFEM_DIR = $(HOME)/PETSC/petscfem-beta-1.82
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test
SRCS = 
NAME = mmove

.PHONY: data force
test_case_data: 
#	-if [ ! "$(cont)" = "y" ] ; then rm -f *.tmp ; fi 
	$(MAKE) data
	octave -q mkmesh.m > /dev/null

data: $(NAME).depl force
$(NAME).depl: force

NPP := $(shell procsel $(PETSCFEM_PROCTABLE))
NP := 1 
test_case: test_case_data
	if [ -e save.state ] ; then cp save.state save.state.tmp ; fi
	-$(MPIRUN) -np $(NP) -machinefile machi.dat $(PROG) \
			-case $(NAME).depl > output.case_$(case).tmp

local_clean::
	-rm -rf mat.output system.dat save.state data.m

LOC_PROG := ns
CASE_NAME := $(NAME)
PROG := $(NS)

#s
