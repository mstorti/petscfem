#__INSERT_LICENSE__
# $Id: Makefile,v 1.1 2001/05/23 20:23:45 mstorti Exp $ 

.PHONY: force all run tests data

all: run

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include ../Makefile.base.test
SRCS = 

#s Test/Advdif tests

NP := $(shell procsel $(PETSCFEM_PROCTABLE))
genload.depl: force
run: genload.depl
	-$(MPIRUN) -np 1 -machinefile machi.dat $(ADVDIF) \
			-case genload.depl #> output.tmp

local_clean::

auxi:
	if [ ! -L advdif ] ; \
	 then ln -sf $(PETSCFEM_DIR)/applications/advdif/advdif.bin advdif ; \
	fi
	echo file advdif >.gdbinit
	echo set args -case genload.depl >>.gdbinit
#s
