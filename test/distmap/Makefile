#__INSERT_LICENSE__
# $Id: Makefile,v 1.26.56.1 2006/02/24 12:26:57 rodrigop Exp $ 

.PHONY: distcont distmap distmat dist_graph

#s Local Targets
#w all tests
tests: distcont distmap distmat dist_graph tests_sb

#w Tests the distributed container class
distcont: distcont.bin
	$(MPIRUN) -np 2     distcont.bin \
		10 100 1e-12   > distcont.sal.tmp

dm2: distmat.bin
	$(MPIRUN) -np 2     distmat.bin \
		4 4 1e-12 0

#w check basic template (map<int,double>)
distmap: distmap.bin
	$(MPIRUN) -np 3     distmap.bin \
		100 10000 1e-12 g > distmapg.sal.tmp
	$(MPIRUN) -np 3     distmap.bin \
		100 10000 1e-12 s > distmaps.sal.tmp

#w check DistMatrix class
distmat: distmat.bin
	$(MPIRUN) -np 3     distmat.bin \
		100 10000 1e-12 0 > distmat.sal.tmp

SRCS = $(wildcard *.cpp)

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
	$(MPIRUN) -np 2     tryme7.bin 20 > dist_graph.np2.out.tmp
	$(MPIRUN) -np 3     tryme7.bin 20 > dist_graph.np3.out.tmp

#w
testsb.bin: testsb.o
	$(MAKE) -C ../../src compile
	-rm -f $@ ; \
	${CXX_CLINKER} -o $@ testsb.o $(LDFLAGS_NO_PROG_LIB)

machi:
	echo minerva >  
	echo minerva >>  
	echo minerva >>  
	echo minerva >>  
	echo minerva >>  

NPSB := 3
#w Tests for `SyncBuffer' and `KeyedOutputBuffer' classes
tests_sb: test_sync_buffer test_kob_plain test_kob_no_sort \
		test_kob_no_key_print test_kob_no_newline

#w Test SyncBuffer class
test_sync_buffer: testsb.bin
	$(MPIRUN) -np $(NPSB) testsb.bin -N 10 -k 0 > $@.output.tmp

#w Test `KeyedOutputBuffer' class
test_kob_plain: testsb.bin
	$(MPIRUN) -np $(NPSB) testsb.bin -k 1 -s 1 -p 1 -n 1 > $@.output.tmp

#w Test `KeyedOutputBuffer' class with no sorting
test_kob_no_sort: testsb.bin
	$(MPIRUN) -np $(NPSB) testsb.bin -k 1 -s 0 -p 1 -n 1 > $@.output.tmp

#w Test `KeyedOutputBuffer' class with no key-printing 
test_kob_no_key_print: testsb.bin
	$(MPIRUN) -np $(NPSB) testsb.bin -k 1 -s 1 -p 0 -n 1 > $@.output.tmp

#w Test `KeyedOutputBuffer' class with no new-lines 
test_kob_no_newline: testsb.bin
	$(MPIRUN) -np $(NPSB) testsb.bin -k 1 -s 1 -p 1 -n 0 > $@.output.tmp

#s
