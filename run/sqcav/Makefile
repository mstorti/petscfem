# $Id: Makefile,v 1.2 2001/01/04 22:45:15 mstorti Exp $ 

.PHONY: all run  compile lclean laplace adv nso laplace_tests force \
        ctest gprof

all: run

gprof: petscfem.prof

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base

NS    = $(PETSCFEM_DIR)/applications/ns/ns.bin

################################################################
NP = 1
PROG  = $(NS)
CASE = sqcav
MCH = -machinefile ./machi.dat

ns: 
	ln -s $(NS) ns

NP := $(shell procsel $(PETSCFEM_PROCTABLE))
run: $(CASE).depl force ns
	$(MAKE) $(PROG)
	-${MPIRUN} -np 1 $(MCH) $(PROG) -case $(CASE).depl

local_clean::
	-chmod +w *.depl
	-rm *.depl *.sal tmp_file*.tmp &> /dev/null
