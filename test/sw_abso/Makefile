default: igrec
#default: test_abso

tests: test_abso

CASE=canal1d

# PETSCFEM_DIR := /u/rodrigop/PETSC/petscfem-hydro
PETSCFEM_DIR := ../..

#NP := $(shell procsel -n $(PETSCFEM_PROCTABLE))
NP := 3

include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

SRCS = 

CLEAN_DIRS = 
DEPEND_DIRS = $(CLEAN_DIRS)

ADVDIF = $(APPDIR)/advdif/advdif_O.bin

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
# Test for abso shallow water 1d

test_abso:
	$(MAKE) $(CASE).depl
	$(OCTAVE) $(CASE).m
	-$(MPIRUN) -np $(NP)  -machinefile ./machi.dat $(ADVDIF) -case $(CASE).depl
	$(OCTAVE)  check_test_abso.m

chiquito:
	$(MAKE) chiquito.depl
	-$(MPIRUN) -np $(NP)  $(ADVDIF) -case chiquito.depl

spiral:
	$(MAKE) spiral.depl
	$(OCTAVE) spiral.m
	-$(MPIRUN) -np $(NP)  $(ADVDIF) -case spiral.depl

corner:
	$(MAKE) corner.depl
	$(OCTAVE) corner.m
	-$(MPIRUN) -np $(NP)  $(ADVDIF) -case corner.depl

igrec:
	$(MAKE) igrec.depl
	$(OCTAVE) igrec.m
	-$(MPIRUN) -np $(NP)  $(ADVDIF) -case igrec.depl


local_clean::
	-chmod +w *.depl
	-rm *.depl *.tmp &> /dev/null

