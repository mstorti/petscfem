# $Id: Makefile,v 1.2 2001/01/24 01:50:20 mstorti Exp $ 

CASE = turbchan
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
$(CASE): $(CASE).depl
	-$(MPIRUN) -np 1 $(ADVDIF) -case $(CASE).depl > $(CASE).out
	verify.pl > turbchan.ver.tmp
