# $Id: Makefile,v 1.2 2001/01/02 16:52:40 mstorti Exp $ 

.PHONY: fstack force

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include ../Makefile.base.test
SRC_DIR = $(PETSCFEM_DIR)/src
PETSCFEM_LIB = $(SRC_DIR)/libpetscfem.a
SRCS = tfstack.cpp

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
# Checks 
fstack: tfstack.bin
	tfstack.bin >tfstack.sal

$(SRC_DIR)/%.o: force
	$(MAKE) -C $(SRC_DIR) $@

$(PETSCFEM_LIB): force
	$(MAKE) -C $(SRC_DIR) compile

tfstack.bin: tfstack.o $(PETSCFEM_LIB)
	${CXX_CLINKER} -o $@ tfstack.o $(SRC_DIR)/fstack.o \
		-lc ${PETSC_SLES_LIB} \
		-L$(LIBRETTO_LIB) -libretto

#
