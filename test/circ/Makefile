#__INSERT_LICENSE__
# $Id: Makefile,v 1.3 2001/04/01 03:11:30 mstorti Exp $ 

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
