#__INSERT_LICENSE__
# $Id: Makefile,v 1.9 2002/07/25 22:35:10 mstorti Exp $ 

#s Local Targets
#w all tests
tests: graph_test

.PHONY: graph_test

tryme: tryme.bin

SRCS = tryme.cpp tryme3.cpp

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

CLEAN_DIRS = 
DEPEND_DIRS = $(CLEAN_DIRS)

LIB = libpetscfem$(osfx).a
LIB_P = ../../src/libpetscfem$(osfx).a

$(LIB_P)(%.o):
	$(MAKE) -C ../../src/ $(LIB)($*.o)

tryme.bin: tryme.o ../../src/libpetscfem.a(graph.o)
	-rm -f $@
	${CXX_CLINKER} -o $@ tryme.o  $(LDFLAGS_NO_PROG_LIB)

# Test Graph class
graph_test: tryme.bin
	sh graph.sh

LDFLAGS = -lglib 

tryme2: tryme2.o
	${CXX_CLINKER} -o $@ tryme2.o  $(LDFLAGS)

tryme3.bin: tryme3.o
	$(MAKE) -C ../../src compile
	${CXX_CLINKER} -o $@ $<  $(LDFLAGS)

run3:
	procsel $(PROCTABLE)
	mpirun -np 2 -machinefile machi.dat tryme3.bin

#s
