# $Id: Makefile,v 1.1 2003/02/17 14:53:21 mstorti Exp $ 

tests: idmap

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
#w Idmap class
idmap: tidmap.bin
	tidmap.bin >tidmap.sal

tidmap.bin: tidmap.o $(LIBPETSCFEM)
	${CXX_CLINKER} -o $@ tidmap.o $(LIBPETSCFEM) \
		-lc ${PETSC_SLES_LIB} -L$(NEWMAT) -lnewmat  $(MESCHACH)/libmes.a

SRCS = 
PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

#s
