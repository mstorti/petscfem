#__INSERT_LICENSE__
# $Id: Makefile,v 1.5 2003/01/08 15:54:26 mstorti Exp $ 

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
	-$(MPIRUN) -np 1 $(ADV) -case circ.dat >/dev/null
	$(OCTAVE) check_circ.m  >check_circ.sal

#
