# $Id: Makefile,v 1.6 2001/01/15 11:11:13 mstorti Exp $ 

.PHONY: all run  compile lclean laplace adv nso laplace_tests force \
        ctest gprof

all: run

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include ../../test/Makefile.base.test

NS = $(PETSCFEM_DIR)/applications/ns/ns.bin
NSO = $(HOME)/PETSC/petscfem-2.8s/applications/ns/ns.bin
#NSO = /home/mstorti/PETSC/petscfem-2.8s/applications/ns/ns.bin 
#NSO = nso.bin

################################################################
NP = 1
ifeq ($(prog),)
$(warning prog not defined)
prog := NS
endif

PROG := $($(prog))
CASE = sqcav
MCH = -machinefile ./machi.dat

NP := $(shell procsel $(PETSCFEM_PROCTABLE))
run: $(CASE).depl force
	@echo starting run
#	$(MAKE) $(PROG)
	-${MPIRUN} -np 1 $(MCH) $(PROG) -case $(CASE).depl

EPERLFLAGS = -d prog=$(prog)
$(CASE).depl: force

local_clean::
	-chmod +w *.depl
	-rm *.depl *.sal tmp_file*.tmp &> /dev/null
