#__INSERT_LICENSE__
# $Id: Makefile,v 1.11 2002/12/12 03:10:08 mstorti Exp $ 

.PHONY: force all run tests data

run: test_case

OLD_MMOVE := 0
export OLD_MMOVE
ifeq ($(OLD_MMOVE),1)
PETSCFEM_DIR = $(HOME)/PETSC/petscfem-mesh-move
else
PETSCFEM_DIR = ../..
endif
# PETSCFEM_DIR = $(HOME)/PETSC/petscfem-beta-1.82
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test
SRCS = 

.PHONY: data force
test_case_data: force;
	$(MAKE) data

ifeq ($(NAME),)
NAME = step3d
endif

LOC_PROG := ns
CASE_NAME := $(NAME)
ifeq ($(OLD_MMOVE),1)
PROG := $(PETSCFEM_DIR)/applications/ns/ns.bin
else
PROG := $(PETSCFEM_DIR)/applications/ns/ns_g.bin
endif

# ifeq ($(NAME),step)
# PROG := $(PETSCFEM_DIR)/applications/ns/ns_g.bin
# endif
step3d:
	$(MAKE) NAME=step3d OCT_SCRIPT=mkstep3d.m test_case

step:
	$(MAKE) NAME=step OCT_SCRIPT=mkstep.m test_case

mmove:
	$(MAKE) NAME=mmove test_case

data: $(NAME).depl force
$(NAME).depl: force

NPP := $(shell procsel $(PETSCFEM_PROCTABLE))
NP := 1 
test_case: test_case_data
	if [ ! "$(OCT_SCRIPT)" = "" ] ; then octave -qH $(OCT_SCRIPT) ; fi
	-$(MPIRUN) -np $(NP) -machinefile machi.dat $(PROG) \
			-case $(NAME).depl

local_clean::
	-rm -rf mat.output system.dat save.state data.m \
		lastmesh.dat remeshing.dat


#s
