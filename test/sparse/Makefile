#__INSERT_LICENSE__
# $Id: Makefile,v 1.10 2002/07/19 20:16:22 mstorti Exp $ 

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

tryme3.bin: tryme3.o
	${CXX_CLINKER} -o $@ $< 

tryme4.bin: tryme4.cpp
	g++ -g -o $@ $<

tryme5.o: CPPFLAGS = -I/u/mstorti/SOFT/boost_1_28_0/ -ftemplate-depth-30
tryme5.bin: tryme5.o
	g++ -g -o $@ $<

#s
