#__INSERT_LICENSE__
# $Id: Makefile,v 1.1 2002/10/06 20:24:15 mstorti Exp $ 
#s Local Targets
#w all tests
tests: pfmat pfmat2

.PHONY: pfmat pfmat2 force iisdmat_petsc iisdmat_slu petscmat \
			spd_petsc spd_slu

LIB_OBJS := 

SRCS = leak.cpp 

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

LIB_OBJS_LIB := $(patsubst %,$(LIB)(%.o),$(LIB_OBJS))

$(LIB)(%.o):
	$(MAKE) -C $(SRC_DIR) 'libpetscfem$(osfx).a($*.o)'

CLEAN_DIRS = 
DEPEND_DIRS = $(CLEAN_DIRS)

#w Build executable
leak.bin: $(LIB_OBJS_LIB) leak.o
	-rm -f $@
	${CXX_CLINKER} -o $@ leak.o  $(LDFLAGS_NO_PROG_LIB)


#s
