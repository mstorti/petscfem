#__INSERT_LICENSE__
# $Id: Makefile,v 1.10 2005/07/05 21:47:01 mstorti Exp $ 

default: swturb

CASE = turbchanw
tests: heat swturb

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

SRCS = 

CLEAN_DIRS = 
DEPEND_DIRS = $(CLEAN_DIRS)

# ADVDIF = $(APPDIR)/advdif/advdif.bin
#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
# Test for shallow water

.PHONY: $(CASE)

data: 
	$(OCTAVE) tcini.m
	$(MAKE) $(CASE).depl

$(CASE).depl: force
$(CASE): data
	-rm -rf some.tmp
	-$(MPIRUN) -np 1 $(ADVDIF) -case $(CASE).depl > $(SUBCASE).sal.tmp
	-mv -f state.out.tmp $(SUBCASE).out.tmp

heat:
	$(MAKE) SUBCASE=heat EPERLFLAGS='-d heat=1' $(CASE)

swturb:
	$(MAKE) SUBCASE=swturb $(CASE)
	$(OCTAVE) check.swturb.m > check.swturb.verif.tmp
