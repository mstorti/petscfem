#__INSERT_LICENSE__
# $Id: Makefile,v 1.15 2001/08/07 16:59:43 mstorti Exp $ 

.PHONY: force all run tests data

run:
	$(MAKE) NP=5 case=iisd case_size=small test_case
	mv -f save.state save.state.iisd.np2.tmp

#w All tests
tests: test_iisd test_partitioning

test_eff: 
	$(MAKE) NP=2 Re=100 case=iisd case_size=large test_case
	mv -f save.state save.state.iisd_peri.np2.tmp

#w Tests IISD with periodic boundary conditions
test_distmap: 
	$(MAKE) NP=2 Re=100 peri=1 case=iisd case_size=small2 test_case
	mv -f save.state save.state.iisd_peri.np2.tmp
	$(MAKE) NP=1 Re=100 peri=1 case=iisd case_size=small2 test_case
	mv -f save.state save.state.iisd_peri.np1.tmp
	octave -q check_peri.m > check_peri.verif.tmp

#w Tests IISD 
test_iisd:
	$(MAKE) NP=2 case=iisd case_size=small test_case
	mv -f save.state save.state.iisd.np2.tmp
	$(MAKE) NP=2 case=iisd ksp=cgs case_size=small test_case
	mv -f save.state save.state.iisd.cgs.np2.tmp
	$(MAKE) NP=1 case=iisd case_size=small test_case
	mv -f save.state save.state.iisd.np1.tmp
	$(MAKE) NP=1 case=lu case_size=small test_case
	mv -f save.state save.state.lu.tmp
	octave -q check_iisd.m > check_iisd.verif.tmp

#w Tests different methods of partitioning. 
test_partitioning: 
	$(MAKE) NP=2 case=metis_part test_case
	mv -f save.state save.state.metis.tmp
	$(MAKE) NP=2  case=hitchhiking_part test_case
	mv -f save.state save.state.hitchhiking.tmp
	$(MAKE) NP=1  case=nearest_neighbor_part test_case
	mv -f save.state save.state.nearest_neighbor.tmp
	octave -q check_part.m > check_part.verif.tmp

PETSCFEM_DIR = ../..
#PETSCFEM_DIR = $(HOME)/PETSC/petscfem-beta-1.82
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test
SRCS = 
NAME = wallke

PROG = $(NS)
LOC_PROG = ns

.PHONY: data force
test_case_data: 
#	-if [ ! "$(cont)" = "y" ] ; then rm -f *.tmp ; fi 
	$(MAKE) data
	octave -q mkmesh.m > /dev/null

data: $(NAME).depl force
$(NAME).depl: force

NPP := $(shell procsel $(PETSCFEM_PROCTABLE))
NPP := 1
ifeq ($(NP),)
NP := $(NPP)
endif
test_case: test_case_data
	if [ -e save.state ] ; then cp save.state save.state.tmp ; fi
	-$(MPIRUN) -np $(NP) -machinefile machi.dat $(PROG) \
			-case $(NAME).depl  > output.case_$(case).tmp

local_clean::
	-rm -rf mat.output system.dat save.state data.m

PROG = $(NS)
LOC_PROG = ns
CASE_NAME = wallke
#s
