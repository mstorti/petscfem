# $Id: Makefile,v 1.1 2001/11/30 00:46:29 mstorti Exp $ 

.PHONY: all run  compile lclean laplace adv nso laplace_tests force \
        ctest gprof

tests: test_re1000

all: run

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

NS = $(PETSCFEM_DIR)/applications/ns/ns.bin

################################################################
NP = 1

PROG := $(NS)
CASE = sqcav
MCH = -machinefile ./machi.dat

NP := $(shell procsel $(PETSCFEM_PROCTABLE))
#w Runs program
#run: $(CASE).depl force
run: data
	-${MPIRUN} -np 1 $(MCH) $(PROG) -case $(CASE).depl  > output.out

test_re1000:
	$(MAKE) run
	rm -f u.Re1000.tmp 
	mv sqcav.state.tmp u.Re1000.tmp 
	octave -q < check_re1000.m > check.re1000.verif.tmp

#w Makes local link for running with debugger
link:
	ln -sf $(PROG) ns

data: $(CASE).depl
	octave -q mksqcav.m

$(CASE).depl: force

local_clean::
	-chmod +w *.depl
	-rm *.depl *.sal tmp_file*.tmp &> /dev/null
