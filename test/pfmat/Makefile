#__INSERT_LICENSE__
# $Id: Makefile,v 1.1.2.1 2001/12/09 14:00:20 mstorti Exp $ 

#s Local Targets
#w all tests
tests: sparse

pfmat: pfmat.bin

SRCS = pfmat.cpp

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

CLEAN_DIRS = 
DEPEND_DIRS = $(CLEAN_DIRS)

SRC_DIR = $(PETSCFEM_DIR)/src
LIB := $(SRC_DIR)/libpetscfem.a
LIB_OBJS := graph util2
LIB_OBJS_LIB := $(patsubst %,$(LIB)(%.o),$(LIB_OBJS))

$(LIB)(%.o):
	$(MAKE) -C $(SRC_DIR) 'libpetscfem.a($*.o)'

pfmat.bin: $(LIB_OBJS_LIB) pfmat.o
	-rm -f $@
	${CXX_CLINKER} -o $@ pfmat.o  $(LDFLAGS) # $(SUPERLU)/SRC/util.o

pfmat: pfmat.bin
	./pfmat.bin
#s
