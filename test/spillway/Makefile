#__INSERT_LICENSE__
# $Id: Makefile,v 1.6 2003/03/23 20:52:43 mstorti Exp $ 

.PHONY: force all run tests data viscous inviscid sphere spillway_mmv

default: spillway

tests: 

spillway:  alehook.efn fifos
	$(MAKE) CASE_NAME=spillway run

spillway_mesh:
	initia=0 octave -qH spillway.m

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

mesh_move:
	$(MAKE) CASE_NAME=spillway_mmv run > mesh_move.output.tmp

alehook.efn:  ../plate/fifo.o
alehook.efn: EFN_LDFLAGS = ../plate/fifo.o

../plate/fifo.o: force
	$(MAKE) -C ../plate fifo.o

.PHONY: force
NP := $(shell procsel $(PETSCFEM_PROCTABLE))
NP := 1
# NP := $(shell procsel ./proctable.tmp)
run: $(CASE_NAME).depl alehook.efn
#	../../tools/eperl_min -P $(CASE_NAME).epl > $(CASE_NAME).depl
	-$(MPIRUN) -np $(NP) -machinefile machi.dat $(PROG) -case $< 

$(CASE_NAME).depl: force

local_clean::
	-rm -rf mat.output *.fifo

#s
