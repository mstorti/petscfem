#__INSERT_LICENSE__
# $Id: Makefile,v 1.4 2007/03/13 01:31:21 mstorti Exp $ 

.PHONY: force all run tests data

all: run

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include ../Makefile.base.test
SRCS = 

#s Test/Advdif tests

NP := $(shell procsel $(PETSCFEM_PROCTABLE))
newff.depl: force
run: newff.depl
	-$(MPIRUN) -machinefile machi.dat -np 2 $(ADVDIF) \
			-case newff.depl </dev/null > output.tmp

local_clean::

LOC_PROG := advdif
PROG := $(ADVDIF)
CASE_NAME := newff
auxi_o:
	if [ ! -L advdif ] ; \
	 then ln -sf $(PETSCFEM_DIR)/applications/advdif/advdif.bin advdif ; \
	fi
	echo file advdif >.gdbinit
	echo set args -case newff.depl >>.gdbinit
#s
