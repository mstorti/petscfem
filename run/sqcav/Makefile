# $Id: Makefile,v 1.3 2001/01/05 23:26:39 mstorti Exp $ 

.PHONY: all run  compile lclean laplace adv nso laplace_tests force \
        ctest gprof

all: run

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base

NS    = $(PETSCFEM_DIR)/applications/ns/ns.bin

################################################################
NP = 1
PROG  = $(NS)
CASE = sqcav
MCH = -machinefile ./machi.dat

ns: $(NS)
	ln -s $(NS) ns

$(NS):
	$(MAKE) -C $(PETSCFEM_DIR) ns

NP := $(shell procsel $(PETSCFEM_PROCTABLE))
run: $(CASE).depl force ns
	echo starting run
	$(MAKE) $(PROG)
	-${MPIRUN} -np 1 $(MCH) $(PROG) -case $(CASE).depl

local_clean::
	-chmod +w *.depl
	-rm *.depl *.sal tmp_file*.tmp &> /dev/null

