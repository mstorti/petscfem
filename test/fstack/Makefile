#__INSERT_LICENSE__
# $Id: Makefile,v 1.5 2001/05/26 19:58:17 mstorti Exp $ 

.PHONY: fstack force

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

../../src/libpetscfem.a(fstack.o): force
	$(MAKE) -C ../../src/ 'libpetscfem.a(fstack.o)'

tfstack.bin: tfstack.o ../../src/libpetscfem.a(fstack.o)
	${CXX_CLINKER} -o $@ tfstack.o $(LIBPETSCFEM) \
		-lc ${PETSC_SLES_LIB} \
		-L$(LIBRETTO_LIB) -libretto

#w Prints all main targets for 'make' along with a brief description
what:
	@$(PETSCFEM_DIR)/tools/makewhat.pl -s Makefile $(OMKFLS)
#
