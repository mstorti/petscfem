#__INSERT_LICENSE__
# $Id: Makefile,v 1.12 2001/08/16 03:54:24 mstorti Exp $ 

#s Local Targets
#w all tests
tests: distmap distmat

dm: distmap.bin
	time $(MPIRUN) -np 6 -machinefile machi.dat distmap.bin \
		1000 10000000 1e-12 g  > distmap.sal.tmp

# check basic template (map<int,double>)
distmap: distmap.bin
	$(MPIRUN) -np 3 -machinefile machi.dat distmap.bin \
		100 10000 1e-12 g > distmapg.sal.tmp
	$(MPIRUN) -np 3 -machinefile machi.dat distmap.bin \
		100 10000 1e-12 s > distmaps.sal.tmp

#w check DistMatrix class
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

#LDFLAGS := -L/usr/local/lib -lc ${PETSC_SLES_LIB}

#w Builds binary 
distmap.bin: distmap.o # libpetscfem
	-rm -f $@
	${CXX_CLINKER} -o $@ distmap.o   $(LDFLAGS)

../../src/libpetscfem.a(%.o): 
	$(MAKE) -C ../../src/ 'libpetscfem.a(%.o)'

#w Builds binary 
distmat.bin: distmat.o ../../src/libpetscfem.a(utils.o)
	-rm -f $@
	${CXX_CLINKER} -o $@ distmat.o $(LIBPETSCFEM) $(LDFLAGS) 

P = tryme
$(P).bin: $(P).cpp
	g++ -o $(P).bin $(P).cpp

#s
