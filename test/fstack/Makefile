# $Id: Makefile,v 1.3 2001/01/20 10:15:29 mstorti Exp $ 

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

tfstack.bin: tfstack.o libpetscfem
	${CXX_CLINKER} -o $@ tfstack.o $(LIBPETSCFEM) \
		-lc ${PETSC_SLES_LIB} \
		-L$(LIBRETTO_LIB) -libretto

#
