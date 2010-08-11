
PETSCFEM_DIR = ../..

PFEM_UP_DIR = ..

include $(PETSCFEM_DIR)/Makefile.base

ADVDIF  = $(PETSCFEM_DIR)/applications/advdif/advdif_O.bin

NAME := plano

eccho:
	echo $(NAME)
	echo $(MPIEXEC)

CASE := $(NAME)
MKMESH := mkplano.m

PROG := $(ADVDIF)
LOC_PROG := advdif

NP := $(shell procsel -n ./proctable)
#NP := 3

default: run


mesh: data
	octave -qH mkplano.m

data: 
	make-depl $(CASE).epl $(CASE).depl 

run: mesh
	mkstepsdir.pl
	make-depl $(CASE).epl $(CASE).depl 
	$(MPIEXEC) -n 1 $(PROG) -case $(CASE).depl
#	time $(MPIEXEC) -1 -machinefile ./machi.dat -n $(NP) $(PROG) -case $(CASE).depl

mpdstart:
	let NP1=$(NP)+1 ; \
	mpdboot --file=./machi.dat -n $$NP1 --rsh=ssh

#w generate DX files
dxo:
	eperl -P viewer.dx > viewer_tmp.dx

dx_case := comp_corner_supersonic_dx
dx:
	$(EPERL) -P -d dx=1 $(dx_case).epl > $(dx_case).depl
	$(MPIRUN) -np 1 $(PETSCFEM_DIR)/applications/advdif/advdif_O.bin \
		 -case $(dx_case).depl

dx_make_command:
	echo hi

rm_all_step:
	rm -rf STEPS* time_history.tmp furn_*

local_clean::
	-chmod +w *.depl
	-rm -rf *.depl *.sal *.tmp *.res STEPS*


