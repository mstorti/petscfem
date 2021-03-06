#__INSERT_LICENSE__

.PHONY: force all run tests data

default:
	$(MAKE) NP=2 case=iisd case_size=small part=random test_case

test_eff: 
	$(MAKE) NP=2 Re=1000 case=iisd_eff case_size=large2 test_case

#w All tests
tests: test_iisd test_partitioning test_distmap \
		test_direct_petsc test_iisd_subpart test_superlu

tests_suite: test_distmap test_direct_petsc 

#w IISD with subpartitioning and partitioning with coalescence
test_iisd_subpart: 
	-if [ -f save.state ] ; then rm -f save.state ; fi
	$(MAKE) NP=1 Re=100 case=lu case_size=small test_case
	mv -f save.state save.state.lu.tmp
	-if [ -f save.state ] ; then rm -f save.state ; fi
	$(MAKE) NP=2 Re=100 sub_part= case=iisd case_size=small test_case
	mv -f save.state save.state.sub_part.tmp
	$(OCTAVE) check_subpart.m > check_subpart.verif.tmp

#w Tests IISD with periodic boundary conditions
test_distmap: 
	$(MAKE) NP=2 Re=100 peri=1 case=iisd case_size=small2 test_case
	mv -f save.state save.state.iisd_peri.np2.tmp
	$(MAKE) NP=1 Re=100 peri=1 case=iisd case_size=small2 test_case
	mv -f save.state save.state.iisd_peri.np1.tmp
	$(OCTAVE) check_peri.m > check_peri.verif.tmp

#w Tests IISD 
test_iisd:
	$(MAKE) NP=2 case=iisd case_size=small part=random test_case
	mv -f save.state save.state.iisd.np2.rand.tmp
	$(MAKE) NP=2 case=iisd case_size=small test_case
	mv -f save.state save.state.iisd.np2.tmp
	$(MAKE) NP=2 case=iisd ksp=cgs case_size=small test_case
	mv -f save.state save.state.iisd.cgs.np2.tmp
	$(MAKE) NP=1 case=iisd case_size=small test_case
	mv -f save.state save.state.iisd.np1.tmp
	$(MAKE) NP=1 case=lu case_size=small test_case
	mv -f save.state save.state.lu.tmp
	$(OCTAVE) check_iisd.m > check_iisd.verif.tmp

#w Tests IISD wirh local_solver = SuperLU
test_iisd_superlu:
ifneq ($(SUPERLU),)
	$(MAKE) NP=2 local_solver=iisd_superlu \
			case=iisd case_size=small part=random test_case
	mv -f save.state save.state.iisd_superlu.np2.rand.tmp
	$(MAKE) NP=2 local_solver=iisd_superlu case=iisd case_size=small test_case
	mv -f save.state save.state.iisd_superlu.np2.tmp
	$(MAKE) NP=2 local_solver=iisd_superlu case=iisd \
			ksp=cgs case_size=small test_case
	mv -f save.state save.state.iisd_superlu.cgs.np2.tmp
	$(MAKE) NP=1 local_solver=iisd_superlu case=iisd case_size=small test_case
	mv -f save.state save.state.iisd_superlu.np1.tmp
	$(MAKE) NP=1 case=lu case_size=small test_case
	mv -f save.state save.state.lu.tmp
endif
	$(OCTAVE) check_iisd_superlu.m > check_iisd_superlu.verif.tmp

#w Tests IISD wirh local_solver = SuperLU
test_superlu:
ifneq ($(SUPERLU),)
	$(MAKE) NP=1 local_solver=direct_superlu case=iisd case_size=small test_case
	mv -f save.state save.state.direct_superlu.tmp
	$(MAKE) NP=1 case=lu case_size=small test_case
	mv -f save.state save.state.lu.tmp
# else
#	$(MAKE) NP=1 local_solver=direct_superlu case=iisd case_size=small test_case
endif
	$(OCTAVE) check_direct_superlu.m > check_direct_superlu.verif.tmp

#w Tests IISD wirh local_solver = SuperLU
test_direct_petsc:
	$(MAKE) NP=1 local_solver=direct_petsc case=iisd case_size=small test_case
	mv -f save.state save.state.direct_petsc.tmp
	$(MAKE) NP=1 case=lu case_size=small test_case
	mv -f save.state save.state.lu.tmp
	$(OCTAVE) check_direct_petsc.m > check_direct_petsc.verif.tmp

#w Tests different methods of partitioning. 
test_partitioning: 
	$(MAKE) NP=2 case=metis_part test_case
	mv -f save.state save.state.metis.tmp
# 	$(MAKE) NP=2  case=hitchhiking_part test_case
# 	mv -f save.state save.state.hitchhiking.tmp
# 	$(MAKE) NP=2  case=nearest_neighbor_part test_case
# 	mv -f save.state save.state.nearest_neighbor.tmp
	$(MAKE) NP=2  case=random_part test_case
	mv -f save.state save.state.random.tmp
	$(OCTAVE) check_part.m > check_part.verif.tmp

PETSCFEM_DIR = ../..
#PETSCFEM_DIR = $(HOME)/PETSC/petscfem-beta-1.82
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test
SRCS = 
NAME = wallke

.PHONY: data force
test_case_data: 
#	-if [ ! "$(cont)" = "y" ] ; then rm -f *.tmp ; fi 
	$(MAKE) data
	$(OCTAVE) mkmesh.m > /dev/null

data: $(NAME).depl force
$(NAME).depl: force

NPP := $(shell procsel $(PETSCFEM_PROCTABLE))
NP := 1 
test_case: test_case_data
	if [ -e save.state ] ;				\
		then cp save.state save.state.tmp ; fi
	-$(MPIEXEC) -n $(NP) $(PROG) -case $(NAME).depl	\
		</dev/null > output.case_$(case).tmp

local_clean::
	-rm -rf mat.output system.dat save.state data.m

LOC_PROG := ns
CASE_NAME := wallke
ifeq ($(PROG),)
PROG := $(NS)
endif

link:
	$(MAKE) -C $(PETSCFEM_DIR) ns
	ln -f $(PROG) ns

#s
