#__INSERT_LICENSE__
# $Id: Makefile,v 1.7 2002/07/18 02:57:33 mstorti Exp $ 

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
	sp.bin PETSc > output.sparse.petsc.tmp
	sp.bin SuperLU > output.sparse.superlu.tmp

tryme.bin: tryme.o
	g++ -g -o $@ $<

tryme2.bin: tryme2.o
	g++ -g -o $@ $<
#s
