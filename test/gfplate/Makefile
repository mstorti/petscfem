# $Id: Makefile,v 1.64 2005/03/23 01:47:51 mstorti Exp $ 

.PHONY: force all run tests default strip gasflow nozzle  vtubei

default: run_shock3d

BOPT := O_c++
SRCS = 
# PETSCFEM_DIR = /u/nnigro/PETSC/petscfem-bubbly9
PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

ifeq ($(CASE_NAME),) 
CASE_NAME := gfshock3d
endif
LOC_PROG := advdif
ifeq ($(PROG),)
PROG := $(ADVDIF)
endif

NP := 1
NL := 

viewer.dx.tmp: viewer.dx
	make-depl $< $@

gbody.efn: gbody.cpp

run_shock3d:
	mkstepsdir.pl
	$(MAKE) CASE_NAME=gfshock3d run

run_shock:
	mkstepsdir.pl
	$(MAKE) CASE_NAME=gfshock run

run_adre:
	$(MPIRUN) $(NL) -np $(NP) -machinefile		\
	machi.dat $(PROG) -case adre.epl

run2d:
	$(MAKE) CASE_NAME=gfabso2d run

runa:
	$(MAKE) CASE_NAME=gfabso run

run_movshock: 
	mkstepsdir.pl
	$(MAKE) CASE_NAME=gfmovshock run

run_abso2dn: gbody.efn
#	mkstepsdir.pl
	$(MAKE) CASE_NAME=gfabso2dn run

run_peri:
	mkstepsdir.pl
	$(MAKE) CASE_NAME=gfperi run

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

mesh: 
	make-depl $(CASE_NAME).epl $(CASE_NAME).depl -D mkmesh=1

run:
	make-depl $(CASE_NAME).epl $(CASE_NAME).depl
	$(MPIRUN) $(NL) -np $(NP) -machinefile			\
	machi.dat $(PROG) -case $(CASE_NAME).depl 

STEPS_DIR := ./STEPS
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
	echo in target dx_make_command: dx_step: $$dx_step
endif
#	state_file=./$(CASE_NAME).dx-state.tmp 
#		octave -qH fixstate.m ; 

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
