
.PHONY:  force data test_plain test_var_eta test_transient rung run runs

SRCS = SRCS = $(wildcard *.cpp)

PETSCFEM_DIR = ../..
BOPT = O_c++
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

#w All tests
tests:  fastmat2a fastmat2b fastmat2bctx2 fastmat2c fastvec

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>------
#w FastMat2 test a 
fastmat2a: testfm2a.bin
	testfm2a.bin > fastmat2a.out.tmp

testfm2a.bin: testfm2a.o $(LIBPETSCFEM)(fastmat2.o 		\
					fm2eperl.o util2.o) 
	g++ $(PROF_FLAGS) -o testfm2a.bin testfm2a.o $(LDFLAGS_NO_PROG_LIB)

testfm2d.bin: testfm2d.o
	g++ $(PROF_FLAGS) -o testfm2d.bin testfm2d.o $(LDFLAGS_NO_PROG_LIB)

effi:
	$(MAKE) testfm2d.bin
	for j in $(seq 1 10) ; \
	do ( ./testfm2d.bin -n 30 -N 5000 -M 1 -t 1 | grep rate ) ; \
	done

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>--
#w FastMat2 test b but with CacheCtx2
fastmat2b: testfm2.bin
	testfm2.bin > fastmat2.out.tmp

testfm2.bin: testfm2.o $(LIBPETSCFEM)(fastmat2.o fm2eperl.o \
		fm2eperl2.o util2.o dgefa.o) 
	${CXX_CLINKER} $(PROF_FLAGS) -o testfm2.bin \
		testfm2.o $(LDFLAGS_NO_PROG_LIB)

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>--
#w FastMat2 test b
fastmat2bctx2: testfm2ctx2.bin
	testfm2ctx2.bin > fastmat2ctx2.out.tmp

testfm2ctx2.bin: testfm2ctx2.o $(LIBPETSCFEM)(fastmat2.o fm2eperl.o \
		fm2eperl2.o util2.o dgefa.o) 
	${CXX_CLINKER} $(PROF_FLAGS) -o testfm2ctx2.bin \
		testfm2ctx2.o $(LDFLAGS_NO_PROG_LIB)

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
