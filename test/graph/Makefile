#__INSERT_LICENSE__
# $Id: Makefile,v 1.8 2002/07/24 15:46:11 mstorti Exp $ 

#s Local Targets
#w all tests
tests: graph

tryme: tryme.bin

SRCS = tryme.cpp tryme3.cpp

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

CLEAN_DIRS = 
DEPEND_DIRS = $(CLEAN_DIRS)

../../src/libpetscfem(%.o):
	$(MAKE) -C ../../src/ libpetscfem.a($*.o)

tryme.bin: tryme.o ../../src/libpetscfem.a(graph.o)
	-rm -f $@
	${CXX_CLINKER} -o $@ tryme.o  $(LDFLAGS) # $(SUPERLU)/SRC/util.o

# Test Graph class
graph: tryme.bin
	sh graph.sh

LDFLAGS = -lglib 

tryme2: tryme2.o
	${CXX_CLINKER} -o $@ tryme2.o  $(LDFLAGS)

LDFLAGS =					\
	  -L$(PETSCFEM_DIR)/src			\
          -lpetscfem$(osfx)			\
          -L$(METIS) -L/usr/local/lib		\
          -lmetis  -L$(NEWMAT)  -lnewmat	\
          -lglib -L$(LIBRETTO_LIB) -libretto	\
          -lANN -L$(ANN_HOME)/lib		\
          -lc ${PETSC_SLES_LIB}			\
          $(MESCHACH)/libmes.a			\
          $(SUPERLULIB)

tryme3.bin: tryme3.o
	$(MAKE) -C ../../src compile
	${CXX_CLINKER} -o $@ $<  $(LDFLAGS)

run3:
	procsel $(PROCTABLE)
	mpirun -np 2 -machinefile machi.dat tryme3.bin

#s
