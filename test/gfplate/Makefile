# $Id: Makefile,v 1.71 2005/05/13 21:33:06 mstorti Exp $ 

.PHONY: force all run tests default strip gasflow nozzle  vtubei

default: run_shock

# BOPT := O_c++
SRCS = 
# PETSCFEM_DIR = /u/nnigro/PETSC/petscfem-bubbly9
PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

ifeq ($(CASE_NAME),) 
#CASE_NAME := gftuyere
CASE_NAME := gfabso
endif
LOC_PROG := advdif
ifeq ($(PROG),)
PROG := $(ADVDIF)
endif

NP := 1
NL := 

viewer.dx.tmp: viewer.dx
	make-depl $< $@

gbody.efn: gbody.cpp

run_shock3d:
	mkstepsdir.pl
	$(MAKE) CASE_NAME=gfshock3d run

run_shock:
	mkstepsdir.pl
	$(MAKE) CASE_NAME=gfshock run

run_adre:
	$(MPIEXEC) $(NL) -n $(NP) -machinefile		\
	machi.dat $(PROG) -case adre.epl

run2d:
	$(MAKE) CASE_NAME=gfabso2d run

runa:
	$(MAKE) CASE_NAME=gfabso run

run_movshock: 
	mkstepsdir.pl
	$(MAKE) CASE_NAME=gfmovshock run

run_abso2dn: gbody.efn
#	mkstepsdir.pl
	$(MAKE) CASE_NAME=gfabso2dn run

run_peri:
	mkstepsdir.pl
	$(MAKE) CASE_NAME=gfperi run

run_cyl: 
	mkstepsdir.pl
	$(MAKE) CASE_NAME=cylabso run

run_cc: 
	mkstepsdir.pl
	$(MAKE) CASE_NAME=comp_corner run

runp: 
	mkstepsdir.pl
	$(MAKE) CASE_NAME=gfplate run

run_nozzle: 
	mkstepsdir.pl
	$(MAKE) CASE_NAME=gfnozzle run

run_nozzle2: 
	mkstepsdir.pl
	$(MAKE) CASE_NAME=gfnozzle2 run

mesh: 
	make-depl $(CASE_NAME).epl $(CASE_NAME).depl -D mkmesh=1

run:
	make-depl $(CASE_NAME).epl $(CASE_NAME).depl
	$(MPIEXEC) $(NL) -machinefile machi.dat -n $(NP) \
	$(PROG) -case $(CASE_NAME).depl 

# STEPS_DIR := ./STEPS
STEPS_DIR := /u/nnigro/PETSC/RUN/COMPRESSION_CORNER/Tuyere/gfplate/RUN_CRANK-NICHOLSON/STEPS
dx_make_command:
ifeq (0,1)
	@state=$(STEPS_DIR)/$(CASE_NAME).state_$$dx_step.tmp ;	\
	if [ -f $$state ] ; then				\
		echo "state file: $$state"  ;			\
		cp $$state ./$(CASE_NAME).dx-state.tmp ;	\
	elif [ -f $$state.gz ] ; then				\
		echo "state file: $$state.gz " ;		\
		gunzip -c $$state.gz >				\
		./$(CASE_NAME).dx-state.tmp ;			\
	else echo "ERROR: NO STATE FILE FOUND!!" ;		\
		echo "tried  $$state and  $$state.gz" ; fi ;	\
	echo in target dx_make_command: dx_step: $$dx_step
endif
#	state_file=./$(CASE_NAME).dx-state.tmp 
#		octave -qH fixstate.m ; 
ifeq (1,1)
	if [ -f gfshock3d.dx-state-3d.tmp ] ; then rm gfshock3d.dx-state-3d.tmp ; fi
	rcp geronimo:$(STEPS_DIR)/gfshock3d.state_$$dx_step.tmp.gz gfshock3d.dx-state-3d.tmp.gz
	octave -qH compstate.m
	rm gfshock3d.dx-state-3d.tmp.gz
endif
ifeq (1,0)
#	cp STEPS/gfnozzle2.state_$$dx_step.tmp ./gfnozzle2.dx-state.tmp 
	@echo "making dx-state for dx_step $$dx_step" ; \
	gunzip -c STEPS/gfnozzle2.state_$$dx_step.tmp.gz > ./gfnozzle2.dx-state.tmp 
endif

dx: 
	make-depl $(CASE_NAME).epl $(CASE_NAME)-dx.depl -d dx=1
	-$(MPIEXEC) $(NL) -n $(NP) -machinefile		\
	machi.dat $(PROG) -case $(CASE_NAME)-dx.depl # > $(CASE_NAME).log

yuv_compress:
	yuv=YUV/tuyere.$$dx_step.yuv ;		\
	tiff=tuyere.tiff ;			\
	if [ -f $$tiff ] ; then			\
		if [ -f $$yuv ] ;		\
			then rm -f $$yuv ;	\
		fi ;				\
		convert $$tiff $$yuv ;		\
		gzip -f $$yuv ;			\
		rm -f $$tiff ;			\
	fi

link:
	$(MAKE) -C ../../ $(LOC_PROG)
	ln -f $(PROG) $(LOC_PROG)

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

#s
