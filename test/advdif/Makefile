# $Id: Makefile,v 1.1.1.1 2000/12/28 12:54:43 mstorti Exp $ 

.PHONY: plano

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include ../Makefile.base.test
SRCS = 

# with this it doesn't work
ADVDIF    = $(PETSCFEM_DIR)/applications/advdif/advdif.bin 

data: burgers.depl plano.con plano.nod plano.m
#data: plano.depl plano.con plano.nod

plano.depl: plano.epl plano.m
burgers.depl: burgers.epl plano.m

plano.con plano.nod: plano.m
	octave -q plano.m


include makefile.d

run: data
	-$(MPIRUN) -np 4 -machinefile machi.dat $(ADVDIF) -case burgers.depl 
