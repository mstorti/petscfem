#__INSERT_LICENSE__
# $Id: Makefile,v 1.21 2003/04/17 22:50:39 mstorti Exp $ 

#s Local Targets
#w all tests
tests: sparse dvector test_cloud

sp: sp.bin

SRCS = sp.cpp tryme7.cpp tryme8.cpp

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

CLEAN_DIRS = 
DEPEND_DIRS = $(CLEAN_DIRS)

LIB = libpetscfem$(osfx).a
LIB_P = ../../src/libpetscfem$(osfx).a

$(LIB_P)(%.o):
	$(MAKE) -C ../../src/ $(LIB)($*.o)

sp.bin: sp.o $(LIB_P)(sparse.o) $(LIB_P)(randomg.o) $(LIB_P)(spsolve.o)
	-rm -f $@
	${CXX_CLINKER} -o $@ sp.o  $(LDFLAGS_NO_PROG_LIB)

sparse: sp.bin
	sp.bin PETSc > output.sparse.petsc.tmp
ifneq ($(SUPERLU),)
	sp.bin SuperLU > output.sparse.superlu.tmp
else
	grep -v '^__' ../sparse.test > output.sparse.superlu.tmp
endif

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
tryme6.o: tryme6.cpp graphs.h
	g++ -c $(CXXOPTFLAGS) -o $@ $< -I$(PETSCFEM_DIR)

tryme6.bin: tryme6.o
	g++ $(CXXOPTFLAGS) -o $@ $< 

#---<*>---//---<*>---//---<*>---//---<*>---//---<*>---// 
dvector: tryme8.bin
	$< > $(basename $<).output.tmp

#w Basic check for dvector class
tryme7.o: tryme7.cpp 
	g++ -c $(CXXOPTFLAGS) -o $@ $< -I$(PETSCFEM_DIR)

tryme7.bin: tryme7.o
	g++ $(CXXOPTFLAGS) -o $@ $< 

#w Checks for dvector class
tryme8.o: tryme8.cpp 
	g++ -c $(CXXOPTFLAGS) -o $@ $< -I$(PETSCFEM_DIR)

tryme8.bin: tryme8.o
	g++ $(CXXOPTFLAGS) -o $@ $< 

#w Checks for cloud class
tryme9.o: tryme9.cpp 
	g++ -c $(CXXOPTFLAGS) $(CPPFLAGS) -o $@ $<

tryme9.bin: tryme9.o $(LIB_P)
	g++ $(CXXOPTFLAGS) -o $@ $< $(LDFLAGS_NO_PROG_LIB)

link9:
	$(MAKE) -C $(PETSCFEM_DIR)/src compile
	$(MAKE) tryme9.bin

test_cloud: link9
	tryme9.bin > tryme9.output.tmp
	octave -qH check_cloud.m > verif.cloud.tmp

#s
