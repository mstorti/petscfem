#__INSERT_LICENSE__
# $Id: Makefile,v 1.5 2002/07/18 02:57:33 mstorti Exp $ 

#s Local Targets
#w all tests
tests: graph

tryme: tryme.bin

SRCS = tryme.cpp

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

#s
