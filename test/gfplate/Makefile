# $Id: Makefile,v 1.27 2005/01/28 17:34:33 mstorti Exp $ 

.PHONY: force all run tests default strip gasflow nozzle  vtubei

default: runc

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

runa2d:
	$(MAKE) CASE_NAME=gfabso2d run

runa:
	$(MAKE) CASE_NAME=gfabso run

runa2dn:
	$(MAKE) CASE_NAME=gfabso2dn run

runc: 
	$(MAKE) CASE_NAME=cylabso run

runp: 
	$(MAKE) CASE_NAME=gfplate run

run: 
	make-depl $(CASE_NAME).epl $(CASE_NAME).depl
	-$(MPIRUN) $(NL) -np $(NP) -machinefile		\
	machi.dat $(PROG) -case $(CASE_NAME).depl # > $(CASE_NAME).log

STEPS_DIR := ./STEPS
## STEPS_DIR := $(HOME)/PETSC/gfplate/STEPS
dx_make_command:
	@state=$(STEPS_DIR)/cylabso.state_$$dx_step.tmp ;	\
	if [ -f $$state ] ; then				\
		echo "state file: $$state"  ;			\
		cp $$state ./cylabso.dx-state.tmp ;		\
	else if [ -f $$state.gz ] ; then			\
		echo "state file: $$state.gz " ;		\
		gunzip -c $$state.gz >				\
		./cylabso.dx-state.tmp ; fi ; fi ;		\
	echo in target dx_make_command: dx_step: $$dx_step

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

#s
