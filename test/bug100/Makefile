#__INSERT_LICENSE__
# $Id: Makefile,v 1.1 2001/05/12 22:33:27 mstorti Exp $ 

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

auxi:
	if [ ! -L advdif ] ; \
	 then ln -sf $(PETSCFEM_DIR)/applications/advdif/advdif.bin advdif ; \
	fi
	echo file advdif >.gdbinit
	echo set args -case newff.depl >>.gdbinit
#s
