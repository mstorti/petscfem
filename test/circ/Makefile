#__INSERT_LICENSE__

.PHONY: constraint_bug

default: tests

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include ../Makefile.base.test
SRCS = 

tests: constraint_bug

# with this it doesn't work
ADVO    = $(PETSCFEM_DIR)/../petscfem-1.48g/applications/advective/adv.bin 

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
# Checks 
constraint_bug: 
	-$(MPIEXEC) -n 1 $(ADV) -case circ.dat $(call in_backg,$@)
	$(OCTAVE) check_circ.m  >check_circ.sal

#
