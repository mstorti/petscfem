# $Id: Makefile,v 1.8 2001/01/19 21:45:52 mstorti Exp $ 

.PHONY: all run  compile lclean laplace adv nso laplace_tests force \
        ctest gprof

all: run

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include ../../test/Makefile.base.test

NS = $(PETSCFEM_DIR)/applications/ns/ns.bin
NSO = $(PETSCFEM_DIR)/../petscfem-2.34g/applications/ns/ns.bin

################################################################
NP = 1

PROG := $(NSO)
CASE = sqcav
MCH = -machinefile ./machi.dat

NP := $(shell procsel $(PETSCFEM_PROCTABLE))
#w Runs program
run: $(CASE).depl force
	-${MPIRUN} -np 1 $(MCH) $(PROG) -case $(CASE).depl

#w Makes local link for running with debugger
link:
	ln -sf $(PROG) ns

data: $(CASE).depl

$(CASE).depl: force

local_clean::
	-chmod +w *.depl
	-rm *.depl *.sal tmp_file*.tmp &> /dev/null
