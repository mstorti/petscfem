tests: test_abso

CASE=canal1d

PETSCFEM_DIR := ../..

NP := $(shell procsel -n $(PETSCFEM_PROCTABLE))
NP := 1

include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

SRCS = 

CLEAN_DIRS = 
DEPEND_DIRS = $(CLEAN_DIRS)

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
# Test for abso shallow water 1d

test_abso:
	$(MAKE) $(CASE).depl
	$(OCTAVE) $(CASE).m
	-$(MPIRUN) -np $(NP)  $(ADVDIF) -case $(CASE).depl > $(CASE).log.tmp
	$(OCTAVE) check_test_abso.m > test_abso.verif.tmp


local_clean::
	-chmod +w *.depl
	-rm *.depl *.sal tmp_file*.tmp &> /dev/null
