#__INSERT_LICENSE__
# $Id: Makefile,v 1.2 2003/02/08 22:15:18 mstorti Exp $ 

#s Local Targets
#w all tests
tests: graph

.PHONY: tryme
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
