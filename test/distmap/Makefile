#__INSERT_LICENSE__
# $Id: Makefile,v 1.4 2001/08/02 01:49:51 mstorti Exp $ 

#s Local Targets
#w all tests
tests: distmap distmat

#w Runs test
distmap: distmap.bin
	$(MPIRUN) -np 3 -machinefile machi.dat distmap.bin \
		100 10000 1e-12 > distmap.sal.tmp

#w Runs test
distmat: distmat.bin
	$(MPIRUN) -np 3 -machinefile machi.dat distmat.bin \
		100 10000 1e-12 0 > distmat.sal.tmp

SRCS = distmap.cpp distmat.cpp

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

CLEAN_DIRS = 
DEPEND_DIRS = $(CLEAN_DIRS)

NPP := $(shell procsel $(PETSCFEM_PROCTABLE))
#w Builds binary 
distmap.bin: distmap.o # libpetscfem
	-rm -f $@
	${CXX_CLINKER} -o $@ distmap.o $(LIBPETSCFEM) $(LDFLAGS)

#w Builds binary 
distmat.bin: distmat.o # libpetscfem
	-rm -f $@
	${CXX_CLINKER} -o $@ distmat.o $(LIBPETSCFEM) $(LDFLAGS)

#s
