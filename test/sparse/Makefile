#__INSERT_LICENSE__
# $Id: Makefile,v 1.5 2001/09/24 20:14:58 mstorti Exp $ 

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
		../../src/libpetscfem.a(randomg.o) \
		../../src/libpetscfem.a(spsolve.o)
	-rm -f $@
	${CXX_CLINKER} -o $@ sp.o  $(LDFLAGS) # $(SUPERLU)/SRC/util.o

sparse: sp.bin
	sp.bin > output.sparse.tmp

tryme.bin: tryme.o
	g++ -g -o $@ $<
#s
