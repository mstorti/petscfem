#__INSERT_LICENSE__

.PHONY: force all run tests data

tests: run
all: run

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include ../Makefile.base.test
SRCS = 

#s Test/Advdif tests
# MCH := -machinefile machi.dat 
newff.depl: force
run: newff.depl
	-$(MPIEXEC) $(MCH) -n 2 $(ADVDIF) \
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
