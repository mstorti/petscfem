#__INSERT_LICENSE__
# $Id: Makefile,v 1.2 2001/09/21 15:37:57 mstorti Exp $ 

#s Local Targets
#w all tests
tests: sparse

sp: sp.bin

SRCS = sp.cpp

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

CLEAN_DIRS = 
DEPEND_DIRS = $(CLEAN_DIRS)

../../src/libpetscfem(%.o):
	$(MAKE) -C ../../src/ libpetscfem.a($*.o)

sp.bin: sp.o ../../src/libpetscfem.a(sparse.o)
	-rm -f $@
	${CXX_CLINKER} -o $@ sp.o   $(LDFLAGS)

tryme.bin: tryme.o
	g++ -g -o $@ $<
#s
