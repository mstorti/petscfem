# $Id: Makefile,v 1.43 2005/02/07 16:04:33 mstorti Exp $ 

.PHONY: force all run tests default strip gasflow nozzle  vtubei

default: run_abso2dn

BOPT := O_c++
SRCS = 
# PETSCFEM_DIR = $(PETSC_ROOT_DIR)/petscfem-beta-3.40
PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

ifeq ($(CASE_NAME),) 
CASE_NAME := cylabso
endif
LOC_PROG := advdif
ifeq ($(PROG),)
PROG := $(ADVDIF)
NL :=
NP := 1
endif

gbody.efn: gbody.cpp

run2d:
	$(MAKE) CASE_NAME=gfabso2d run

runa:
	$(MAKE) CASE_NAME=gfabso run

run_abso2dn: gbody.efn
	$(MAKE) CASE_NAME=gfabso2dn run

run_cyl: 
	mkstepsdir.pl
	$(MAKE) CASE_NAME=cylabso run

run_cc: 
	mkstepsdir.pl
	$(MAKE) CASE_NAME=comp_corner run

runp: 
	mkstepsdir.pl
	$(MAKE) CASE_NAME=gfplate run

run_nozzle: 
	mkstepsdir.pl
	$(MAKE) CASE_NAME=gfnozzle run

run: 
	make-depl $(CASE_NAME).epl $(CASE_NAME).depl
	-$(MPIRUN) $(NL) -np $(NP) -machinefile		\
	machi.dat $(PROG) -case $(CASE_NAME).depl # > $(CASE_NAME).log

## STEPS_DIR := ./STEPS
STEPS_DIR := ./STEPS-consistent-mass
## STEPS_DIR := $(HOME)/PETSC/gfplate/STEPS
dx_make_command:
ifeq (1,1)
	@state=$(STEPS_DIR)/$(CASE_NAME).state_$$dx_step.tmp ;	\
	if [ -f $$state ] ; then				\
		echo "state file: $$state"  ;			\
		cp $$state ./$(CASE_NAME).dx-state.tmp ;	\
	elif [ -f $$state.gz ] ; then				\
		echo "state file: $$state.gz " ;		\
		gunzip -c $$state.gz >				\
		./$(CASE_NAME).dx-state.tmp ;			\
	else echo "ERROR: NO STATE FILE FOUND!!" ;		\
		echo "tried  $$state and  $$state.gz" ; fi ;	\
	state_file=./$(CASE_NAME).dx-state.tmp			\
		octave -qH fixstate.m ;				\
	echo in target dx_make_command: dx_step: $$dx_step
endif

dx: 
	make-depl $(CASE_NAME).epl $(CASE_NAME)-dx.depl -d dx=1
	-$(MPIRUN) $(NL) -np $(NP) -machinefile		\
	machi.dat $(PROG) -case $(CASE_NAME)-dx.depl # > $(CASE_NAME).log

link:
	$(MAKE) -C ../../ $(LOC_PROG)
	ln -f $(PROG) $(LOC_PROG)

zip:
	for f in STEPS*/*.state*.tmp ;		\
		do echo gzipping $$f ;		\
		gzip $$f ;			\
	done

clean_steps:
	rrm -rf STEPS*

#s
