# $Id: Makefile,v 1.4 2001/12/02 15:55:46 mstorti Exp $ 

.PHONY: all run  compile lclean laplace adv nso laplace_tests force \
        ctest gprof

tests: test_re1000

all: run

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

NS = $(PETSCFEM_DIR)/applications/ns/ns.bin

LOC_PROG := ns
CASE_NAME := sqcav
PROG := $(NS)

################################################################
NP = 1

PROG := $(NS)
CASE = sqcav
MCH = -machinefile ./machi.dat

NP := $(shell procsel $(PETSCFEM_PROCTABLE))
#w Runs program
#run: $(CASE).depl force
run:  data
	-${MPIRUN} -np 1 $(MCH) $(PROG) -case $(CASE).depl  > output.out

# tests the square cavity in a 20x20 mesh at Re=1000
# with respect with a reference solution obtained also with Petscfem and
# verified with Ghia results
#w Square cavity test in 20x20 mesh
test_re1000:
	$(MAKE) case=weak_form_0 run
	$(MAKE) case=weak_form_1 run
	octave -q < check_re1000.m > check.re1000.verif.tmp

#w Tests IISD solver on the square cavity problem
test_iisd:
	$(MAKE) case=lu run
	$(MAKE) case=iisd_sbp1 run
	$(MAKE) case=iisd_sbp2 run
	$(MAKE) case=iisd_sbp4 run
	$(MAKE) case=iisd_sbp8 run
	$(MAKE) case=iisd_sbp16 run
	octave -q < check_iisd.m > check.iisd.verif.tmp

#w Makes local link for running with debugger
link:
	ln -sf $(PROG) ns

data: $(CASE).depl
	octave -q mksqcav.m

# EPERLFLAGS = -d case=$$case
$(CASE).depl: force

local_clean::
	-chmod +w *.depl
	-rm *.depl *.sal tmp_file*.tmp &> /dev/null
