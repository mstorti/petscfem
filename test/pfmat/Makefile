#__INSERT_LICENSE__
# $Id: Makefile,v 1.1.2.2 2001/12/24 04:00:44 mstorti Exp $ 

#s Local Targets
#w all tests
tests: pfmat

SRCS = pfmat.cpp

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

CLEAN_DIRS = 
DEPEND_DIRS = $(CLEAN_DIRS)

SRC_DIR = $(PETSCFEM_DIR)/src
LIB := $(SRC_DIR)/libpetscfem.a
LIB_OBJS := graph util2 iisdmat sparse spsolve distmat iisdcr \
		iisdgraph spline seval debug
LIB_OBJS_LIB := $(patsubst %,$(LIB)(%.o),$(LIB_OBJS))

$(LIB)(%.o):
	$(MAKE) -C $(SRC_DIR) 'libpetscfem.a($*.o)'

pfmat.bin: $(LIB_OBJS_LIB) pfmat.o
	-rm -f $@
	${CXX_CLINKER} -o $@ pfmat.o  $(LDFLAGS)

NP := $(shell procsel $(PETSCFEM_PROCTABLE))

pfmat: # pfmat.bin
	$(MPIRUN) -np 2 -machinefile machi.dat ./pfmat.bin

#s
