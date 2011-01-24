 #$Id: Makefile,v 1.5 2007/01/28 00:29:19 mstorti Exp $

.PHONY: force all run run_struct run_mmv

default: tests

tests: test_bdf_dgcl_gaussian_diag test_bdf_dgcl

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

NP := 1
MCH := 

test_bdf_dgcl_gaussian_diag:
	make-depl $(CASE).epl $(CASE).depl -d adv_case=gaussian_diag
	$(OCTAVE) mkmesh.m
	$(MAKE) run $(call redir,$@.output.tmp)
	$(MAKE) mmv-force.efn
	mkstepsdir.pl
	$(MPI_HOME1)/bin/mpiexec $(MCH) -a $(CASE) \
		-np $(NP) $(PROG) -code bdf -case $(CASE).depl \
		$(call redir,$@.output.tmp)
	$(OCTAVE) verif.m > $@.check.tmp


test_bdf_dgcl:
	make-depl $(CASE).epl $(CASE).depl -d adv_case=bdf_dgcl
	$(OCTAVE) mkmesh.m
	$(MAKE) run $(call redir,$@.output.tmp)
	$(MAKE) mmv-force.efn
	mkstepsdir.pl
	$(MPI_HOME1)/bin/mpiexec $(MCH) -a $(CASE) \
		-np $(NP) $(PROG) -code bdf -case $(CASE).depl 
		$(call redir,$@.output.tmp)
	$(OCTAVE) verif.m > $@.check.tmp

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

local_clean::
	-shopt -s nullglob ; rm -rf STEPS* *.vtk
