# $Id: Makefile,v 1.1 2001/02/09 23:54:31 mstorti Exp $ 

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
#	octave -q mkmesh.m

$(CASE).depl: data $(CASE).epl

$(CASE): $(CASE).depl
	-$(MPIRUN) -np 1 $(ADVDIF) -case $(CASE).depl
