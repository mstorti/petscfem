#__INSERT_LICENSE__
# $Id: Makefile,v 1.6 2002/11/03 01:31:22 mstorti Exp $ 

# I have lost `thash2.cpp'!!
tests: texthash # texthash_find

#s Local Targets
#w Runs the test
texthash: tthash.bin
	$< >tthash.sal

texthash_find:: thash2.bin
	$< > thash2.output.tmp

SRCS = tthash.cpp
LIB_OBJS := texthash fstack

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

CLEAN_DIRS = 
DEPEND_DIRS = $(CLEAN_DIRS)

#w Builds the binary of the test
tthash.bin: $(LIB_OBJS_LIB)  tthash.o 
	-rm -f $@
	${CXX_CLINKER} -o $@ tthash.o $(LIBPETSCFEM) $(LDFLAGS_NO_PROG_LIB)

thash2.bin: $(LIB_OBJS_LIB)  thash2.o 
	-rm -f $@
	${CXX_CLINKER} -o $@ thash2.o $(LIBPETSCFEM) $(LDFLAGS_NO_PROG_LIB)


#s
