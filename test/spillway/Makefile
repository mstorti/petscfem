#__INSERT_LICENSE__
# $Id: Makefile,v 1.13 2003/04/02 18:05:19 mstorti Exp $ 

.PHONY: force all run tests data viscous inviscid sphere spillway_mmv

default: wave

spillway:  alehook.efn fifos
	$(MAKE) CASE_NAME=spillway run

spillway_mmv:
	$(MAKE) CASE_NAME=spillway_mmv run > spillway_mmv.output.tmp

spillway_mesh:
	initia=0 octave -qH spillway.m

copy:
	cp wave_mmv.state.tmp.orig  wave_mmv.state.tmp
	cp wave.state.tmp.orig wave.state.tmp

wavei:
	$(MAKE) initia=1 wave

wave:  alehook.efn fifos
	$(MAKE) CASE_NAME=wave run

wave_mmv:
	$(MAKE) CASE_NAME=wave_mmv run > wave_mmv.output.tmp

SRCS = inviscid.cpp
PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

LOC_PROG := ns
PROG := $(NS)

fifos:
	if [ -e ns2mmv.fifo ] ; then rm ns2mmv.fifo ; fi
	mkfifo ns2mmv.fifo
	if [ -e mmv2ns.fifo ] ; then rm mmv2ns.fifo ; fi
	mkfifo mmv2ns.fifo

alehook.efn:  ../plate/fifo.o
alehook.efn: EFN_LDFLAGS = ../plate/fifo.o $(LIBPETSCFEM)

../plate/fifo.o: force
	$(MAKE) -C ../plate fifo.o

.PHONY: force
NP := $(shell procsel -l $(PETSCFEM_PROCTABLE))
NP := 1
run: $(CASE_NAME).depl alehook.efn
#	../../tools/eperl_min -P $(CASE_NAME).epl > $(CASE_NAME).depl
	-$(MPIRUN) -np $(NP) -machinefile machi.dat $(PROG) -case $< 

$(CASE_NAME).depl: force

local_clean::
	-rm -rf mat.output *.fifo

#s
