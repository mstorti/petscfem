#__INSERT_LICENSE__
# $Id: Makefile,v 1.22 2004/01/09 17:03:16 mstorti Exp $ 

.PHONY: distcont distmap distmat dist_graph

#s Local Targets
#w all tests
tests: distcont distmap distmat dist_graph

#w Tests the distributed container class
distcont: distcont.bin
	$(MPIRUN) -np 2 -machinefile machi.dat distcont.bin \
		10 100 1e-12   > distcont.sal.tmp

dm2: distmat.bin
	$(MPIRUN) -np 2 -machinefile machi.dat distmat.bin \
		4 4 1e-12 0

#w check basic template (map<int,double>)
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

#LDFLAGS := -L/usr/local/lib -lc ${PETSCFEM_PETSC_LIBS}

#w Builds binary 
distmap.bin: distmap.o 
	-rm -f $@
	${CXX_CLINKER} -o $@ distmap.o   $(LDFLAGS_NO_PROG_LIB)

#w Builds binary
distcont.bin: distcont.o 
	-rm -f $@
	${CXX_CLINKER} -o $@ distcont.o   $(LDFLAGS_NO_PROG_LIB)

distcont.o: ../../src/distcont.h

LIB := $(SRC_DIR)/libpetscfem$(osfx).a
LIB_OBJS_LIB := $(patsubst %,$(LIB)(%.o),$(LIB_OBJS))

$(LIB)(%.o):
	$(MAKE) -C $(SRC_DIR) 'libpetscfem$(osfx).a($*.o)'

#w Builds binary 
distmat.bin: distmat.o $(LIB)(utils.o) $(LIB)(distmat.o)
	-rm -f $@
	${CXX_CLINKER} -o $@ distmat.o $(LIBPETSCFEM) $(LDFLAGS_NO_PROG_LIB) 

P := tryme4
$(P).bin: $(P).o $(LIBPETSCFEM) 
	-rm $@
	${CXX_CLINKER} $(PROF_FLAGS) -o $@ $(P).o $(OTHER_OBJS) $(LDFLAGS)

tryme5.bin: tryme5.cpp
	g++ -g -o tryme5.bin tryme5.cpp

tryme6.bin: tryme6.cpp
	g++ -g -o tryme6.bin tryme6.cpp

SRC_DIR = $(PETSCFEM_DIR)/src
LIB := $(SRC_DIR)/libpetscfem$(osfx).a
LIB_OBJS := iisdgraph graph debug
LIB_OBJS_LIB := $(patsubst %,$(LIB)(%.o),$(LIB_OBJS))

$(LIB)(%.o):
	$(MAKE) -C $(SRC_DIR) libpetscfem$(osfx).a($*.o)

tryme7.bin: tryme7.o $(LIB_OBJS_LIB)
	-rm -f $@
	${CXX_CLINKER} -o $@ tryme7.o $(LDFLAGS_NO_PROG_LIB)

#w Checks distributed graph class 'StoreGraph'
dist_graph: tryme7.bin
	tryme7.bin 20 > dist_graph.np1.out.tmp
	$(MPIRUN) -np 2 -machinefile machi.dat tryme7.bin 20 > dist_graph.np2.out.tmp
	$(MPIRUN) -np 3 -machinefile machi.dat tryme7.bin 20 > dist_graph.np3.out.tmp

#w
syncbuff.bin: syncbuff.o
	-rm -f $@ ; \
	${CXX_CLINKER} -o $@ syncbuff.o $(LDFLAGS_NO_PROG_LIB)

machi:
	echo minerva > machi.dat
	echo minerva >> machi.dat
	echo minerva >> machi.dat
	echo minerva >> machi.dat
	echo minerva >> machi.dat

#s
