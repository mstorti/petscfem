#__INSERT_LICENSE__
# $Id: Makefile,v 1.2 2001/07/28 20:04:04 mstorti Exp $ 

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
	-$(MPIRUN) -np 2 -machinefile machi.dat $(ADVDIF) \
			-case newff.depl > output.tmp

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
