#__INSERT_LICENSE__
# $Id: Makefile,v 1.13 2002/07/21 19:01:48 mstorti Exp $ 

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
tryme4.o: tryme4.cpp dvector.h
	g++ -c $(CXXOPTFLAGS) -o $@ $< -I$(PETSCFEM_DIR)

tryme4.bin: tryme4.o
	g++ $(CXXOPTFLAGS) -o $@ $< -L$(LIBRETTO_LIB) -libretto 

tryme5.o: CPPFLAGS = -I/u/mstorti/SOFT/boost_1_28_0/ -ftemplate-depth-30
tryme5.bin: tryme5.o
	g++ -g -o $@ $<

CXXOPTFLAGS = -g
tryme6.o: tryme6.cpp dvector.h graphdv.h graphs.h
	g++ -c $(CXXOPTFLAGS) -o $@ $< -I$(PETSCFEM_DIR)

tryme6.bin: tryme6.o
	g++ $(CXXOPTFLAGS) -o $@ $< 

#s
