# $Id: Makefile,v 1.16 2005/01/25 00:57:09 mstorti Exp $ 

.PHONY: force all run tests default strip gasflow nozzle  vtubei

default: runc

BOPT := O_c++
SRCS = 
PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

ifeq ($(CASE_NAME),) 
# CASE_NAME := gfplate
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

runc: 
	$(MAKE) CASE_NAME=cylabso run

run: 
	make-depl $(CASE_NAME).epl $(CASE_NAME).depl
	-$(MPIRUN) $(NL) -np $(NP) -machinefile		\
	machi.dat $(PROG) -case $(CASE_NAME).depl # > $(CASE_NAME).log

runa: 
	make-depl $(CASE_NAME).epl $(CASE_NAME).depl
	echo -n > $(CASE_NAME).some-rslt.tmp
	-$(MPIRUN) $(NL) -np $(NP) -machinefile		\
	machi.dat $(PROG) -case $(CASE_NAME).depl # > $(CASE_NAME).log

dx_make_command:
	cp STEPS/cylabso.state_$$dx_step.tmp ./cylabso.dx-state.tmp
	@echo dx_step: $$dx_step

dx: 
	make-depl $(CASE_NAME).epl $(CASE_NAME)-dx.depl -d dx=1
	-$(MPIRUN) $(NL) -np $(NP) -machinefile		\
	machi.dat $(PROG) -case $(CASE_NAME)-dx.depl # > $(CASE_NAME).log

link:
	$(MAKE) -C ../../ $(LOC_PROG)
	ln -f $(PROG) $(LOC_PROG)

#s
