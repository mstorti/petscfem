# $Id: Makefile,v 1.10 2007/03/13 01:31:21 mstorti Exp $ 

.PHONY:  force data test_plain test_var_eta test_transient rung run runs

SRCS = SRCS = $(wildcard *.cpp)

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

#w All tests
tests:  fastmat2a fastmat2b fastmat2c fastvec

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>------
#w FastMat2 test a 
fastmat2a: testfm2a.bin
	testfm2a.bin > fastmat2a.out.tmp

testfm2a.bin: testfm2a.o $(LIBPETSCFEM)(fastmat2.o 		\
					fm2eperl.o util2.o) 
	g++ $(PROF_FLAGS) -o testfm2a.bin testfm2a.o $(LDFLAGS_NO_PROG_LIB)

testfm2d.bin: testfm2d.o
	g++ $(PROF_FLAGS) -o testfm2d.bin testfm2d.o $(LDFLAGS_NO_PROG_LIB)

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>--
#w FastMat2 test b
fastmat2b: testfm2.bin
	testfm2.bin > fastmat2.out.tmp

testfm2.bin: testfm2.o $(LIBPETSCFEM)(fastmat2.o fm2eperl.o \
		fm2eperl2.o util2.o dgefa.o) 
	${CXX_CLINKER} $(PROF_FLAGS) -o testfm2.bin \
		testfm2.o $(LDFLAGS_NO_PROG_LIB)

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>--
#w FastMat2 test c 
fastmat2c: testfm2c.bin
	testfm2c.bin > fastmat2c.out.tmp

testfm2c.bin: testfm2c.o $(LIBPETSCFEM)(fastmat2.o fm2eperl.o util2.o) 
	g++ $(PROF_FLAGS) -o testfm2c.bin testfm2c.o $(LDFLAGS_NO_PROG_LIB)

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
#w FAST Library
fastvec: tfastvec.bin
	tfastvec.bin > tfastvec.out.tmp

tfastvec.bin: tfastvec.o
	g++ -o tfastvec.bin tfastvec.o 

tfastvec.o: tfastvec.cpp

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
#w Eigenvalues
fm2eig.bin: fm2eig.o 
	g++ -o $@ fm2eig.o $(LDFLAGS_NO_PROG) -lblas -llapack -lg2c 
