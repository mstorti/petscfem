# $Id: Makefile,v 1.6 2005/01/22 22:10:20 mstorti Exp $ 

.PHONY: force all run tests default strip gasflow nozzle  vtubei

default: runa

# BOPT := O_c++
SRCS = 
PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

ifeq ($(CASE_NAME),) 
# CASE_NAME := gfplate
CASE_NAME := gfabso
endif
LOC_PROG := advdif
ifeq ($(PROG),)
PROG := $(ADVDIF)
NL :=
NP := 1
endif

run: 
	make-depl $(CASE_NAME).epl $(CASE_NAME).depl
	-$(MPIRUN) $(NL) -np $(NP) -machinefile		\
	machi.dat $(PROG) -case $(CASE_NAME).depl # > $(CASE_NAME).log

runa: 
	make-depl $(CASE_NAME).epl $(CASE_NAME).depl
	echo -n > gfabso.some-rslt.tmp
	-$(MPIRUN) $(NL) -np $(NP) -machinefile		\
	machi.dat $(PROG) -case $(CASE_NAME).depl # > $(CASE_NAME).log

dx: 
	make-depl $(CASE_NAME).epl $(CASE_NAME)-dx.depl -d dx=1
	-$(MPIRUN) $(NL) -np $(NP) -machinefile		\
	machi.dat $(PROG) -case $(CASE_NAME)-dx.depl # > $(CASE_NAME).log

link:
	$(MAKE) -C ../../ $(LOC_PROG)
	ln -f $(PROG) $(LOC_PROG)

#s
