# $Id: Makefile,v 1.11 2006/08/07 20:36:58 mstorti Exp $ 

.PHONY: force all run tests default strip gasflow nozzle  vtubei

default: run

tests: test_alebcconv

BOPT := g_c++
SRCS = 

PETSCFEM_DIR := ../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test


CASE_NAME := gfabso
LOC_PROG := advdif
PROG := $(ADVDIF)

ifeq ($(HOSTNAME),galileo)
NP := 1
MCH := 
endif

run: mvbody.efn	
	@echo "running with ALE= $(ALE)"
	rm -f gfabso.state.tmp
	make-depl $(CASE_NAME).epl $(CASE_NAME).depl -d use_ALE=$(ALE)
	mkstepsdir.pl
	cp $(CASE_NAME).depl $(CASE_NAME).epl ./STEPS
	$(MPI_HOME)/bin/mpiexec $(MCH) -a $(CASE_NAME)		\
		-np $(NP) $(PROG) -case $(CASE_NAME).depl	\
		$(call redir,output.ALE$(ALE).tmp)

ALE := 1
test_alebcconv:
	$(MAKE) ALE=0 run
	mv -f gfabso.state.tmp gfabso.state-ALE0.tmp
	$(MAKE) ALE=1 run
	mv -f gfabso.state.tmp gfabso.state-ALE1.tmp
	$(OCTAVE) checkinv.m > output.alebcconv.verif.tmp

prog:
	$(MAKE) -C $(PETSCFEM_DIR) $(LOC_PROG)

link:
	$(MAKE) -C $(PETSCFEM_DIR) $(LOC_PROG)
	ln -sf $(ADVDIF_G) .

zip:
	@for f in STEPS*/*.state*.tmp ;		\
		do echo gzipping $$f ;		\
		gzip $$f ;			\
	done

clean_steps:
	rrm -rf STEPS*

gf_shell_hook_init:

gf_shell_hook_pre: zip

gf_shell_hook_post: 

gf_shell_hook_close:

yuv_compress:
	echo "Doing nothing"

jpg:
	@shopt -s nullglob ;					\
	for f in *.tiff ; do					\
		echo "Processing $$f" ;  			\
		ff=`basename $$f .tiff` ;			\
		convert $$ff.tiff $$ff.jpg ;			\
		mv -f $$ff.jpg ~/TEX/papers/dynabc-data ;	\
		ad $$ff.tiff ;					\
	done

kill:
	$(MPI_HOME)/bin/mpdkilljob -a $(CASE_NAME)

#s
