# $Id: Makefile,v 1.3 2001/12/02 14:21:37 mstorti Exp $ 

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
	echo "In target run: case = $$case"
	-${MPIRUN} -np 1 $(MCH) $(PROG) -case $(CASE).depl  > output.out

test_re1000:
	$(MAKE) case=weak_form_0 run
	$(MAKE) case=weak_form_1 run
	octave -q < check_re1000.m > check.re1000.verif.tmp

test_iisd:
	$(MAKE) case=lus  run
	rm -f u.Re1000.tmp 
	mv sqcav.state.tmp u.Re1000.wf0.tmp 
#
	$(MAKE) weak_form=1 run
	rm -f u.Re1000.tmp 
	mv sqcav.state.tmp u.Re1000.wf1.tmp 
#
	octave -q < check_re1000.m > check.re1000.verif.tmp

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
