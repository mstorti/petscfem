# $Id: Makefile,v 1.1.1.1 2000/12/28 12:54:43 mstorti Exp $ 

CASE = turbchan
run: $(CASE)

include ../Makefile.base.test

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base

SRCS = 

CLEAN_DIRS = 
DEPEND_DIRS = $(CLEAN_DIRS)

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
# Test for shallow water

.PHONY: $(CASE)
$(CASE): adv $(CASE).depl
	-$(MPIRUN) -np 1 $(ADV) -case $(CASE).depl # $(CASE).out
