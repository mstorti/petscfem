#__INSERT_LICENSE__
# $Id: Makefile,v 1.2 2002/10/06 21:51:35 mstorti Exp $ 
#s Local Targets
#w all tests
tests: pfmat pfmat2

.PHONY: pfmat pfmat2 force iisdmat_petsc iisdmat_slu petscmat \
			spd_petsc spd_slu

LIB_OBJS := 

SRCS = leak.cpp 

PETSCFEM_DIR = ../..
# PETSCFEM_DIR = $(HOME)/PETSC/petscfem-beta-3.03

include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

ifneq ($(PETSC_VERSION_USED_BY_PETSCFEM),2.1.3)
CPPFLAGS += -DUSE_OLD_PETSC_VERSION
endif

LIB_OBJS_LIB := $(patsubst %,$(LIB)(%.o),$(LIB_OBJS))

$(LIB)(%.o):
	$(MAKE) -C $(SRC_DIR) 'libpetscfem$(osfx).a($*.o)'

CLEAN_DIRS = 
DEPEND_DIRS = $(CLEAN_DIRS)
leak.o: 
#w Build executable
leak.bin: $(LIB_OBJS_LIB) leak.o
	-rm -f $@
	${CXX_CLINKER} -o $@ leak.o  $(LDFLAGS_NO_PROG_LIB)


#s
