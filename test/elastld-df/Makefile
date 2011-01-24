#$Id: Makefile,v 1.1 2006/03/11 23:45:43 mstorti Exp $

.PHONY: force all run run_struct

default: run

PETSCFEM_DIR := ../..
# BOPT := O_c++
SRCS =
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test
MPI_HOME := /usr

CASE := elastld
PROG := $(NS)
LOC_PROG := ns

ifeq ($(HOSTNAME),aquiles)
NP := $(shell procsel -n ./proctable)
MCH := -1 -machinefile ./machi.dat
else
NP := 1
MCH := 
endif

tests:
	$(MAKE) test_ld_elast

test_ld_elast:
	$(MAKE) mesh
	$(MAKE) run $(call redir,$@.output.tmp)
	$(OCTAVE) verif.m > test_ld_elast.check.tmp

mesh:
	make-depl $(CASE).epl $(CASE).depl -dmkmesh=1

split:
	$(PETSCFEM_DIR)/tools/hexasplit.bin \
		-i elastld.con.tmp -o elastld.con-tetra.tmp

run: bdry.efn
	mkstepsdir.pl
	make-depl $(CASE).epl $(CASE).depl -dmkmesh=0
	cp -f data.pl *.epl *.depl ./STEPS/
	$(MPIEXEC) -a $(CASE) $(MCH) -np $(NP) $(PROG) \
		-code struct -case $(CASE).depl

dx:
	make-depl $(CASE).epl $(CASE)-dx.depl -d dx=1
	$(PROG) -case $(CASE)-dx.depl 

zip2:
	@shopt -s nullglob ;				\
	for f in STEPS/$(CASE).state-*.tmp		\
		STEPS/$(CASE).mmv-state-*.tmp	\
		STEPS/flutter-plate.state*.tmp ;	\
	    do echo gzipping $$f ;			\
		gzip -f $$f ;				\
	done

zip:
	shopt -s nullglob ;				\
	find STEPS* -name '$(CASE)*.state*.tmp' |	\
		xargs -r gzip -f 


shell_init:

shell_pre:

shell_post: zip

shell_close:

DIRE := ./STEPS
# DIRE := ./SAVED-STEPS/STEPS-2006-FEB-06-COUPLED-KFAC0-1
dx_make_command:
	gunzip -c $(DIRE)/$(CASE).state-$$dx_step.tmp.gz \
		> $(CASE).dx-state.tmp
	gunzip -c $(DIRE)/flutter.mmv-state-$$dx_step.tmp > \
		flutter.coords-dx-state.tmp

yuv_compress:
	yuv=YUV/$(CASE).$$dx_step.yuv ;		\
	tiff=$(CASE).tiff ;			\
	if [ -f $$tiff ] ; then			\
		if [ -f $$yuv ] ;		\
			then rm -f $$yuv ;	\
		fi ;				\
		convert $$tiff $$yuv ;		\
		gzip -f $$yuv ;			\
		rm -f $$tiff ;			\
	fi

link:
	$(MAKE) -C $(PETSCFEM_DIR) ns
	ln -sf $(NS) ns

bdry.efn: bdry.cpp

mpdstart:
	let NP1=$(NP)+1 ; \
	$(MPI_HOME1)/bin/mpdboot --file=./machi.dat -n $$NP1 --rsh=rsh

kill:
	-mpdkilljob -a $(CASE) &> /dev/null

ltag:
	makeltag save-log.txt notes.txt

local_clean::
	shopt -s nullglob ; rm -f *.fifo ns advdif *.log *.vtk
	rrm -rf STEPS*

yuv-clean:
	shopt -s nullglob ; rrm -f YUV/*.png YUV/*.yuv.gz

video:
	mkavi.sh elastld

prog:
	$(MAKE) -C $(PETSCFEM_DIR) ns

RDIRE := /u/mstorti/PETSC/petscfem-cases/elastld/test2-df/STEPSc
cpy:
	rsync --delete -av aquilesp:$(RDIRE) .
