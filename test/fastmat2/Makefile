# $Id: Makefile,v 1.3 2002/11/03 12:07:22 mstorti Exp $ 

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

testfm2a.bin: testfm2a.o $(LIBPETSCFEM)(fastmat2.o fm2eperl.o util2.o) 
	g++ $(PROF_FLAGS) -o testfm2a.bin testfm2a.o $(LIBPETSCFEM) \
		-L$(NEWMAT) -lnewmat ${PETSC_SLES_LIB} $(MESCHACH)/libmes.a

testfm2d.bin: testfm2d.o
	g++ $(PROF_FLAGS) -o testfm2d.bin testfm2d.o $(LIBPETSCFEM) \
		-L$(NEWMAT) -lnewmat ${PETSC_SLES_LIB} $(MESCHACH)/libmes.a -lblas

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>--
#w FastMat2 test b
fastmat2b: testfm2.bin
	testfm2.bin > fastmat2.out.tmp

testfm2.bin: testfm2.o $(LIBPETSCFEM)(fastmat2.o fm2eperl.o util2.o) 
	${CXX_CLINKER} $(PROF_FLAGS) -o testfm2.bin testfm2.o \
		$(LIBPETSCFEM) ${PETSC_SLES_LIB} -L$(NEWMAT) -lnewmat \
		$(MESCHACH)/libmes.a

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>--
#w FastMat2 test c 
fastmat2c: testfm2c.bin
	testfm2c.bin > fastmat2c.out.tmp

testfm2c.bin: testfm2c.o $(LIBPETSCFEM)(fastmat2.o fm2eperl.o util2.o) 
	g++ $(PROF_FLAGS) -o testfm2c.bin testfm2c.o $(LIBPETSCFEM) \
		 -L$(NEWMAT) -lnewmat ${PETSC_SLES_LIB}  $(MESCHACH)/libmes.a

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
#w FAST Library
fastvec: tfastvec.bin
	tfastvec.bin > tfastvec.out.tmp

tfastvec.bin: tfastvec.o
	g++ -o tfastvec.bin tfastvec.o 

tfastvec.o: tfastvec.cpp
