# $Id: Makefile,v 1.2 2001/01/02 16:52:34 mstorti Exp $ 

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
	octave -q check_circ.m  >check_circ.sal

#
