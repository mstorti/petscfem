 #$Id: Makefile,v 1.5 2007/01/28 00:29:19 mstorti Exp $

.PHONY: force all run run_struct run_mmv

default: tests

tests:
	test_bdf_dgcl

## activar para compilar la opcion sin debager
# BOPT := O_c++
PETSCFEM_DIR := ../..
SRCS = 
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

CASE := gascont
CASE_NAME := $(CASE)
CASE_NAME_MMV := $(CASE)-mmv
PROG := $(ADVDIF)
LOC_PROG := advdif

ifeq ($(HOSTNAME),aquiles)
NP := $(shell procsel -n ./proctable)
MCH := -1 -machinefile ./machi.dat
else
NP := 1
MCH := 
endif

test_bdf_dgcl:
	$(MAKE) mesh
	$(MAKE) run

mesh:
	make-depl $(CASE).epl $(CASE).depl -d mkmesh=1
#	$(OCTAVE) mkmesh.m

run:
	$(MAKE) mmv-force.efn
	mkstepsdir.pl
	make-depl $(CASE).epl $(CASE).depl
	cp -f data.pl *.epl *.depl ./STEPS/
	$(MPI_HOME1)/bin/mpiexec $(MCH) -a $(CASE) \
		-np $(NP) $(PROG) -code bdf -case $(CASE).depl

zip:
	@shopt -s nullglob ;				\
	for f in STEPS/$(CASE).state-*.tmp		\
		STEPS/$(CASE_NAME_MMV).state-*.tmp	\
		STEPS/$(CASE)-str.state*.tmp ;		\
	    do echo gzipping $$f ;			\
		gzip -f $$f ;				\
	done

shell_init:

shell_pre:

shell_post: zip
#	@octave -qH checkstr.m

shell_close:

link:
	$(MAKE) prog
	ln -sf $(PROG) .

prog:
	$(MAKE) -C $(PETSCFEM_DIR) $(LOC_PROG)

kill:
	mpdkilljob -a $(CASE)
