# $Id cases-1.0.31-132-ga72f976 Thu Dec 6 16:33:22 2007 -0300$ 

.PHONY: all run  compile lclean laplace adv nso laplace_tests force \
        ctest gprof srfgath

default: run
tests: test_gasflow_1elem_jaco

# BOPT := O_c++
PETSCFEM_DIR := ../../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

LOC_PROG := advdif
CASE := gfnwttest
PROG := $(ADVDIF)

################################################################
## Default
NP := 1
NL := 
MCH := 

mesh:
	make-depl $(CASE).epl $(CASE).depl -d mkmesh=1

#w Runs program, use `mkmesh=0' or 1 to control mesh creation
test_gasflow_1elem_jaco: 
	make-depl $(CASE).epl $(CASE).depl -d mkmesh=1
	$(MPI_HOME)/bin/mpiexec -a $(CASE) $(NL) $(MCH) \
		-np $(NP) $(PROG) -case $(CASE).depl </dev/null > $@_output.tmp
	$(OCTAVE) $@_verif.m &> $@_verif.tmp

prog:
	$(MAKE) -C $(PETSCFEM_DIR) advdif

local_sw::
	chmod 755 mkmpeg*.sh

local_clean::
	-rm -rf STEPS* 
	-rm -f *.depl *.sal tmp_file*.tmp mat.output &> /dev/null

#w Make link for debugging
link:
	$(MAKE) BOPT=g_c++ -C $(PETSCFEM_DIR) advdif
	ln -sf $(ADVDIF_G) .

kill:
	mpdkilljob -a $(CASE)
