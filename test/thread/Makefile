#__INSERT_LICENSE__
# $Id: Makefile,v 1.1 2001/11/28 02:17:31 mstorti Exp $ 

#s Local Targets
#w all tests
tests: graph

tryme: tryme.bin

SRCS = tryme.cpp

include $(PETSCFEM_DIR)/Makefile.base

CLEAN_DIRS = 
DEPEND_DIRS = $(CLEAN_DIRS)

$(PETSCFEM_DIR)/src/libpetscfem(%.o):
	$(MAKE) -C $(PETSCFEM_DIR)/src/ libpetscfem.a($*.o)

LDFLAGS := -lc -lpthread

tryme.bin: tryme.o 
	-rm -f $@
	${CXX_CLINKER} -o $@ tryme.o  $(LDFLAGS) 

#s
