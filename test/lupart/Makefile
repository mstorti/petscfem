#__INSERT_LICENSE__
# $Id: Makefile,v 1.8 2001/07/25 17:27:31 mstorti Exp $ 

.PHONY: force all run tests data

#w All tests
tests: test_iisd test_partitioning

#w Tests 
test_iisd: 
	$(MAKE) NP=2 case=iisd case_size=small test_case_poiseuille
	mv -f save.state save.state.iisd.np2.tmp
	$(MAKE) NP=2 case=iisd ksp=cgs case_size=small test_case_poiseuille
	mv -f save.state save.state.iisd.cgs.np2.tmp
	$(MAKE) NP=1 case=iisd case_size=small test_case_poiseuille
	mv -f save.state save.state.iisd.np1.tmp
	$(MAKE) NP=1 case=lu case_size=small test_case_poiseuille
	mv -f save.state save.state.lu.tmp
	octave -q check_iisd.m > check_iisd.verif.tmp

#w Tests different methods of partitioning. 
test_partitioning: 
	$(MAKE) NP=2  test_case_metis_part
	mv -f save.state save.state.metis.tmp
	$(MAKE) NP=2  test_case_hitchhiking_part
	mv -f save.state save.state.hitchhiking.tmp
	$(MAKE) NP=1  test_case_nearest_neighbor_part
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
test_case_%_data: 
#	-if [ ! "$(cont)" = "y" ] ; then rm -f *.tmp ; fi 
	$(MAKE) EPERLFLAGS='-d case=$* -d cont=$(cont)' data
	octave -q mkmesh.m > /dev/null

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
			-case $(NAME).depl > output.case_$*.tmp
#			-case $(NAME).depl 
#			-case $(NAME).depl | tee output.case_$*.tmp

local_clean::
	-rm -rf mat.output system.dat save.state data.m

auxi:
	if [ ! -L advdif ] ; \
	 then ln -sf $(PROG) $(LOC_PROG) ; \
	fi
	echo file ns >.gdbinit
	echo set args -case $(NAME).depl >>.gdbinit
#s
