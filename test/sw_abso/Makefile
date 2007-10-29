#__INSERT_LICENSE__
# $Id: Makefile,v 1.11 2004/09/25 11:32:33 mstorti Exp $ 

default: igrec
#default: test_abso

tests: test_abso

CASE=canal1d

# PETSCFEM_DIR := /u/rodrigop/PETSC/petscfem-hydro
PETSCFEM_DIR := ../..

NP := $(shell procsel -n $(PETSCFEM_PROCTABLE))
NP := 3

include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

SRCS = 

CLEAN_DIRS = 
DEPEND_DIRS = $(CLEAN_DIRS)

# ADVDIF = $(APPDIR)/advdif/advdif_O.bin

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
# Test for abso shallow water 1d

test_abso:
	$(MAKE) $(CASE).depl
	$(OCTAVE) $(CASE).m
	-$(MPIEXEC) -n $(NP)  -machinefile ./machi.dat $(ADVDIF) -case $(CASE).depl
	$(OCTAVE)  check_test_abso.m

chiquito:
	$(MAKE) chiquito.depl
	-$(MPIEXEC) -n $(NP)  $(ADVDIF) -case chiquito.depl

spiral:
	$(MAKE) spiral.depl
	$(OCTAVE) spiral.m
	-$(MPIEXEC) -n $(NP)  $(ADVDIF) -case spiral.depl

corner:
	$(MAKE) corner.depl
	$(OCTAVE) corner.m
	-$(MPIEXEC) -n $(NP)  $(ADVDIF) -case corner.depl

igrec:
	$(MAKE) igrec.depl
	$(OCTAVE) igrec.m
	-$(MPIEXEC) -n $(NP)  $(ADVDIF) -case igrec.depl


local_clean::
	-chmod +w *.depl
	-rm *.depl *.tmp &> /dev/null

