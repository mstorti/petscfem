#__INSERT_LICENSE__
# $Id: Makefile,v 1.2 2001/07/31 03:01:35 mstorti Exp $ 

#s Local Targets
#w Runs the test
distmap: distmap.bin
	$< >distmap.sal

SRCS = distmap.cpp

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

CLEAN_DIRS = 
DEPEND_DIRS = $(CLEAN_DIRS)

distmap.o: distmap.cpp

#w Builds the binary of the test

NPP := $(shell procsel $(PETSCFEM_PROCTABLE))
distmap.bin: distmap.o # libpetscfem
	-rm -f $@
	${CXX_CLINKER} -o $@ distmap.o $(LIBPETSCFEM) $(LDFLAGS)

#s
