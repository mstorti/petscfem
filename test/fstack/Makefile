#__INSERT_LICENSE__
# $Id: Makefile,v 1.14 2007/02/20 14:57:47 mstorti Exp $ 

.PHONY: fstack force tests test_autostr

default: fstack

tests: fstack test_autostr

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include ../Makefile.base.test
SRC_DIR = $(PETSCFEM_DIR)/src
PETSCFEM_LIB = $(SRC_DIR)/libpetscfem.a
SRCS = tfstack.cpp

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
#w test for fstack class
fstack: tfstack.bin
	tfstack.bin >tfstack.out
	perl fscheck.pl < fstack2.out.tmp > fstack2.verif.tmp

$(SRC_DIR)/%.o: force
	$(MAKE) -C $(SRC_DIR) $@

tfstack.bin: tfstack.o $(LIBPETSCFEM)
	${CXX_CLINKER} -o $@ tfstack.o $(LDFLAGS_NO_PROG_LIB)

test_autostr: tautostr.bin
	$< > test_autostr.output.tmp

tautostr.bin: tautostr.o $(LIB)(autostr.o)
	${CXX_CLINKER}	-g -Wl,-rpath,$(LD_LIBRARY_PATH) \
		-o $@ tautostr.o $(LIBPETSCFEM) -lc 

#
