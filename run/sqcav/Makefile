# $Id: Makefile,v 1.6.2.1 2001/01/19 16:41:50 mstorti Exp $ 

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

ifeq ($(prog),)
$(warning "Not defined 'prog' variable!! Setting to NS ...")
prog = NS
endif

CASE = sqcav
MCH = -machinefile ./machi.dat

NP := $(shell procsel $(PETSCFEM_PROCTABLE))
EPERLFLAGS = -d prog=$(prog)
run: $(CASE).depl force
	-${MPIRUN} -np 1 $(MCH) $($(prog)) -case $(CASE).depl

$(CASE).depl: force

local_clean::
	-chmod +w *.depl
	-rm *.depl *.sal tmp_file*.tmp &> /dev/null
