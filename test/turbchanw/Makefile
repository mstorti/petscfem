# $Id: Makefile,v 1.2 2001/02/10 15:16:31 mstorti Exp $ 

CASE = turbchanw
run: $(CASE)

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

SRCS = 

CLEAN_DIRS = 
DEPEND_DIRS = $(CLEAN_DIRS)

ADVDIF = $(APPDIR)/advdif/advdif.bin
#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
# Test for shallow water

.PHONY: $(CASE)

data: 
	octave -q tcini.m
	$(MAKE) $(CASE).depl

$(CASE): data
	-$(MPIRUN) -np 1 $(ADVDIF) -case $(CASE).depl
