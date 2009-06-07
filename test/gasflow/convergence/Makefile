# $Id cases-1.0.31-132-ga72f976 Thu Dec 6 16:33:22 2007 -0300$ 

.PHONY: all run  compile lclean laplace adv nso laplace_tests force \
        ctest gprof srfgath

default: run

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

ifeq ($(HOSTNAME),aquiles)
NP := $(shell procsel -n ./proctable)
NL := -1
MCH := -machinefile ./machi.dat
NP := 2
endif

ifeq ($(HOSTNAME),coyote)
NP := $(shell procsel -n ./proctable.coyote)
NL := -1
MCH := -machinefile ./machi.dat
endif

mesh:
	make-depl $(CASE).epl $(CASE).depl -d mkmesh=1

#w Runs program, use `mkmesh=0' or 1 to control mesh creation
run: 
	mkstepsdir.pl
	make-depl $(CASE).epl $(CASE).depl -d mkmesh=0
	$(MPI_HOME)/bin/mpiexec -a $(CASE) $(NL) $(MCH) \
		-np $(NP) $(PROG) -case $(CASE).depl </dev/null 
	$(OCTAVE) proc.m

prog:
	$(MAKE) -C $(PETSCFEM_DIR) advdif

local_sw::
	chmod 755 mkmpeg*.sh

local_clean::
	-rm -rf STEPS* 
	-rm -f *.depl *.sal tmp_file*.tmp # &> /dev/null

shell_init:

shell_pre:

shell_post: 
	$(MAKE) zip

shell_close:

zip:
	shopt -s nullglob ;				\
	find STEPS* -name 'gfnwttest*.state*.tmp' |	\
		xargs -r gzip -f 

#w Make link for debugging
link:
	$(MAKE) BOPT=g_c++ -C $(PETSCFEM_DIR) advdif
	ln -sf $(ADVDIF_G) .

kill:
	mpdkilljob -a $(CASE)
