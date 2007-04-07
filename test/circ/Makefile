#__INSERT_LICENSE__

.PHONY: constraint_bug

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include ../Makefile.base.test
SRCS = 

# with this it doesn't work
ADVO    = $(PETSCFEM_DIR)/../petscfem-1.48g/applications/advective/adv.bin 

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
# Checks 
constraint_bug: 
	-$(MPIRUN) -np 1 $(ADV) -case circ.dat $(call in_backg,$@)
	$(OCTAVE) check_circ.m  >check_circ.sal

#
