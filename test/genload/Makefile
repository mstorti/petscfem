#__INSERT_LICENSE__
# $Id: Makefile,v 1.2 2001/05/24 17:40:02 mstorti Exp $ 

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
	-rm -rf mat.output

auxi:
	if [ ! -L advdif ] ; \
	 then ln -sf $(PETSCFEM_DIR)/applications/advdif/advdif.bin advdif ; \
	fi
	echo file advdif >.gdbinit
	echo set args -case genload.depl >>.gdbinit
#s
