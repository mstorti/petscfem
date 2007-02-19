#__INSERT_LICENSE__
# $Id: Makefile,v 1.13.98.1 2007/02/19 03:35:08 mstorti Exp $ 

.PHONY: fstack force tests test_autostr

default: fstack

#tests: fstack test_autostr
tests:  test_autostr

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

tfstack.bin: tfstack.o $(LIB)(fstack.o)
#	${CXX_CLINKER} -shared -Wl,-rpath,$(LD_LIBRARY_PATH) 
	${CXX_CLINKER}	-g -Wl,-rpath,$(LD_LIBRARY_PATH) -o $@ tfstack.o	\
		$(LIBPETSCFEM) -lc ${PETSCFEM_PETSC_LIBS}			\
		-L$(LIBRETTO_LIB) -libretto

test_autostr: tautostr.bin
	$< > test_autostr.output.tmp

tautostr.bin: tautostr.o $(LIB)(autostr.o)
	${CXX_CLINKER}	-g -Wl,-rpath,$(LD_LIBRARY_PATH) -o $@ tautostr.o	\
		$(LIBPETSCFEM) -lc 

#
