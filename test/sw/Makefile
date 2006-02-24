tests: test_flux

CASE=centerflux

# PETSCFEM_DIR := /u/rodrigop/PETSC/petscfem-hydro
PETSCFEM_DIR := ../..

NP := $(shell procsel -n $(PETSCFEM_PROCTABLE))
NP := 1

include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

SRCS = 

CLEAN_DIRS = 
DEPEND_DIRS = $(CLEAN_DIRS)

# ADVDIF = $(APPDIR)/advdif/advdif_g.bin

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
# Test for abso shallow water 1d

test_flux:
	$(MAKE) $(CASE).depl
	$(OCTAVE) $(CASE).m
	-$(MPIRUN) -np $(NP)     $(ADVDIF) -case $(CASE).depl
#	$(OCTAVE)  check_test_abso.m

local_clean::
	-chmod +w *.depl
	-rm *.depl *.tmp &> /dev/null
