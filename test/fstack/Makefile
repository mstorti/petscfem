#__INSERT_LICENSE__
# $Id: Makefile,v 1.7 2003/02/12 19:32:13 mstorti Exp $ 

.PHONY: fstack force tests

tests: fstack

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
		$(LIBPETSCFEM) -lc ${PETSC_SLES_LIB}			\
		-L$(LIBRETTO_LIB) -libretto

#
