#__INSERT_LICENSE__

.PHONY: distcont distmap distmat dist_graph

#s Local Targets
#w all tests
tests: distcont distmap distmat dist_graph tests_sb

# MCH := -machinefile machi.dat

#w Tests the distributed container class
distcont: distcont.bin
	$(MPIEXEC) $(MCH) -n 2 distcont.bin \
		10 100 1e-12  $(call redir,distcont.sal.tmp)

dm2: distmat.bin
	$(MPIEXEC) $(MCH) -n 2 distmat.bin \
		4 4 1e-12 0 $(call redir,$@)

#w check basic template (map<int,double>)
distmap: distmap.bin
	$(MPIEXEC) $(MCH) -n 3 distmap.bin \
		100 10000 1e-12 g $(call redir,distmapg.sal.tmp)
	$(MPIEXEC) $(MCH) -n 3 distmap.bin \
		100 10000 1e-12 s $(call redir,distmaps.sal.tmp)

#w check DistMatrix class
distmat: distmat.bin
	$(MPIEXEC) $(MCH) -n 3  distmat.bin \
		100 10000 1e-12 0 $(call redir,distmat.sal.tmp)

SRCS = $(wildcard *.cpp)

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

CLEAN_DIRS = 
DEPEND_DIRS = $(CLEAN_DIRS)

#LDFLAGS := -L/usr/local/lib -lc ${PETSCFEM_PETSC_LIBS}

## Redirection from /dev/null to avoid error from mpiexec on output
## Set to blank to run interactive
## redir := 

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
	$(MPIEXEC) $(MCH) -n 2 tryme7.bin 20 \
		$(call redir,dist_graph.np2.out.tmp)
	$(MPIEXEC) $(MCH) -n 3 tryme7.bin 20 \
		$(call redir,dist_graph.np3.out.tmp)

#w
testsb.bin: testsb.o
	$(MAKE) -C ../../src compile
	-rm -f $@ ; \
	${CXX_CLINKER} -o $@ testsb.o $(LDFLAGS_NO_PROG_LIB)

machi:
	echo minerva > machi.dat
	echo minerva >> machi.dat
	echo minerva >> machi.dat
	echo minerva >> machi.dat
	echo minerva >> machi.dat

NPSB := 3
#w Tests for `SyncBuffer' and `KeyedOutputBuffer' classes
tests_sb: test_sync_buffer test_kob_plain test_kob_no_sort \
		test_kob_no_key_print test_kob_no_newline

#w Test SyncBuffer class
test_sync_buffer: testsb.bin
	$(MPIEXEC) -n $(NPSB) testsb.bin -Ne 10 -k 0 $(RDR) \
		$(call redir,$@.output.tmp)

#w Test `KeyedOutputBuffer' class
test_kob_plain: testsb.bin
	$(MPIEXEC) -n $(NPSB) testsb.bin -k 1 -s 1 -p 1 -n 1 \
		$(call redir,$@.output.tmp)

#w Test `KeyedOutputBuffer' class with no sorting
test_kob_no_sort: testsb.bin
	$(MPIEXEC) -n $(NPSB) testsb.bin -k 1 -s 0 -p 1 -n 1 \
		$(call redir,$@.output.tmp)

#w Test `KeyedOutputBuffer' class with no key-printing 
test_kob_no_key_print: testsb.bin
	$(MPIEXEC) -n $(NPSB) testsb.bin -k 1 -s 1 -p 0 -n 1 \
		$(call redir,$@.output.tmp)

#w Test `KeyedOutputBuffer' class with no new-lines 
test_kob_no_newline: testsb.bin
	$(MPIEXEC) -n $(NPSB) testsb.bin -k 1 -s 1 -p 1 -n 0 $(RDR) \
		$(call redir,$@.output.tmp)

#s
