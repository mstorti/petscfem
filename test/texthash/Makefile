# $Id: Makefile,v 1.1 2001/03/03 08:21:03 mstorti Exp $ 

#s Local Targets
#w Runs the test
texthash: tthash.bin
	$< >tthash.sal

SRCS = tthash.cpp

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

CLEAN_DIRS = 
DEPEND_DIRS = $(CLEAN_DIRS)

#w Builds the binary of the test
tthash.bin: tthash.o libpetscfem
	-rm -f $@
	${CXX_CLINKER} -o $@ tthash.o $(LIBPETSCFEM) $(LDFLAGS)

#s
