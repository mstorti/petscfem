#__INSERT_LICENSE__
# $Id: Makefile,v 1.3 2001/09/22 04:29:44 mstorti Exp $ 

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

sp.bin: sp.o ../../src/libpetscfem.a(sparse.o) \
		../../src/libpetscfem.a(randomg.o)
	-rm -f $@
	${CXX_CLINKER} -o $@ sp.o   $(LDFLAGS)

tryme.bin: tryme.o
	g++ -g -o $@ $<
#s
