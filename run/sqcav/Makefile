# $Id: Makefile,v 1.1 2001/01/04 20:13:05 mstorti Exp $ 

.PHONY: all run  compile lclean ns laplace adv nso laplace_tests force \
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
run: $(CASE).depl force
	procsel proctable
	$(MAKE) $(PROG)
	-${MPIRUN} -np $(NP) $(MCH) $(PROG) -case $(CASE).depl

local_clean::
	-chmod +w *.depl
	-rm *.depl *.sal tmp_file*.tmp &> /dev/null
