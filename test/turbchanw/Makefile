#__INSERT_LICENSE__

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
	-$(MPIEXEC) -n 1 $(ADVDIF) -case $(CASE).depl \
		</dev/null > $(SUBCASE).sal.tmp
	-mv -f state.out.tmp $(SUBCASE).out.tmp

heat:
	$(MAKE) SUBCASE=heat EPERLFLAGS='-d heat=1' $(CASE)

swturb:
	$(MAKE) SUBCASE=swturb $(CASE)
	$(OCTAVE) check.swturb.m > check.swturb.verif.tmp
