#__INSERT_LICENSE__
# $Id: Makefile,v 1.19 2003/10/11 14:03:44 mstorti Exp $ 

.PHONY: force all run tests data

default: step

tests: step step3d

run: test_case

PETSCFEM_DIR = ../..
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
PROG := $(NS)

#w Deformation of cube with a axisymmetric ramp on one face
step3d:
	@$(MAKE) -C ../../tools hexasplit.bin >/dev/null
	$(MAKE) NAME=step3d OCT_SCRIPT=mkstep3d.m test_case > step3d.out.tmp
	$(OCTAVE) proc2.m > step3d.verif.tmp

#w Deformation of square with a ramp on one side
step:
	$(MAKE) NAME=step OCT_SCRIPT=mkstep.m test_case > step.out.tmp
	$(OCTAVE) proc.m > step.verif.tmp

#w Check whether lockers are created correctly in adaptor
locker:
	$(MAKE) NAME=step case=locker OCT_SCRIPT=mkstep.m \
		test_case > step-locker.out.tmp
	$(OCTAVE) proc.m > step.verif.tmp

mmove:
	$(MAKE) NAME=mmove test_case

data: $(NAME).depl force
$(NAME).depl: force

NPP := $(shell procsel $(PETSCFEM_PROCTABLE))
NP := 1 
test_case: test_case_data
	if [ ! "$(OCT_SCRIPT)" = "" ] ; then $(OCTAVE) $(OCT_SCRIPT) ; fi
	-$(MPIRUN) -np $(NP) -machinefile machi.dat $(PROG) \
			-case $(NAME).depl

local_clean::
	-rm -rf mat.output system.dat save.state data.m \
		lastmesh.dat remeshing.dat


#s
