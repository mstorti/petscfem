tests: test_abso

CASE=canal1d

PETSCFEM_DIR := /u/mstorti/PETSC/petscfem-hydro

NP := $(shell procsel -n $(PETSCFEM_PROCTABLE))
NP := 1

include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

SRCS = 

CLEAN_DIRS = 
DEPEND_DIRS = $(CLEAN_DIRS)

ADVDIF = $(APPDIR)/advdif/advdif_g.bin
#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
# Test for abso shallow water 1d

test_abso:
	$(MAKE) $(CASE).depl
	octave -q $(CASE).m
	-$(MPIRUN) -np $(NP)  $(ADVDIF) -case $(CASE).depl
	$(octave) -qH check_test_abso.m


local_clean::
	-chmod +w *.depl
	-rm *.depl *.sal tmp_file*.tmp &> /dev/null
