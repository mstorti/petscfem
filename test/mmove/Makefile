#__INSERT_LICENSE__
# $Id: Makefile,v 1.14 2002/12/26 14:36:13 mstorti Exp $ 

.PHONY: force all run tests data

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

step3d:
	@$(MAKE) -C ../../tools hexasplit.bin >/dev/null
	$(MAKE) NAME=step3d OCT_SCRIPT=mkstep3d.m test_case > step3d.out.tmp
	octave -qH proc2.m > step3d.verif.tmp

step:
	$(MAKE) NAME=step OCT_SCRIPT=mkstep.m test_case > step.out.tmp
	octave -qH proc.m > step.verif.tmp

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
