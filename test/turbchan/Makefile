#__INSERT_LICENSE__
# $Id: Makefile,v 1.5 2002/07/25 22:35:10 mstorti Exp $ 

CASE = turbchan
#s Local Targets
#w [default] Runs the test
run: $(CASE)

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

SRCS = 

CLEAN_DIRS = 
DEPEND_DIRS = $(CLEAN_DIRS)

# ADVDIF = $(APPDIR)/advdif/advdif.bin
#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
# Test for shallow water

.PHONY: $(CASE)
$(CASE): $(CASE).depl
	-$(MPIEXEC) -n 1 $(ADVDIF) -case $(CASE).depl > $(CASE).out
	verify.pl > turbchan.ver.tmp

#s
