# $Id: Makefile,v 1.9 2001/12/07 18:27:30 mstorti Exp $ 

.PHONY: all run  compile lclean laplace adv nso laplace_tests force \
        ctest gprof

pp:
	$(MAKE) case=iisd_sbp1_uj runs

tests: test_re1000 test_iisd

all: run

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

NS = $(PETSCFEM_DIR)/applications/ns/ns.bin

LOC_PROG := ns
#CASE_NAME := qharm
CASE_NAME := sqcav
PROG := $(NS)

################################################################
NP = 1

PROG := $(NS)
MCH = -machinefile ./machi.dat

NP := $(shell procsel $(PETSCFEM_PROCTABLE))

#w Runs program
run:  data
	-${MPIRUN} -np 1 $(MCH) $(PROG) -case $(CASE).depl | tee  output.out

export CASE
runs: CASE := sqcav
runs:
	$(MAKE) run

runq: CASE := qharm
runq: force
#	$(MAKE) case=iisd_sbp1 run
	for n in 1 2 4 6 8 16 32 64 128 ; \
		do time $(MAKE) case=iisd_sbp$$n run ; \
		mv output.out output.iisd$$n.out ; \
	done

runs: CASE := sqcav
test_big:
	time $(MAKE) case=big runs

# tests the square cavity in a 20x20 mesh at Re=1000
# with respect with a reference solution obtained also with Petscfem and
# verified with Ghia results
#w Square cavity test in 20x20 mesh
test_re1000:
	$(MAKE) case=weak_form_0 runs
	$(MAKE) case=weak_form_1 runs
	octave -q < check_re1000.m > check.re1000.verif.tmp

#w Tests IISD solver on the square cavity problem
test_iisd:
	$(MAKE) case=lu runs
	$(MAKE) case=iisd_sbp1 runs
	$(MAKE) case=iisd_sbp2 runs
	$(MAKE) case=iisd_sbp4 runs
	$(MAKE) case=iisd_sbp8 runs
	$(MAKE) case=iisd_sbp16 runs
	octave -q < check_iisd.m > check.iisd.verif.tmp

#w Makes local link for running with debugger
link:
	ln -sf $(PROG) ns

data: force
	echo in target data: CASE = $(CASE)
	$(MAKE) $(CASE).depl
	echo "in target data : CASE $(CASE)"
	octave -q mk$(CASE).m

# EPERLFLAGS = -d case=$$case
$(CASE).depl: force
#	echo in target \$(CASE).depl : CASE $(CASE)

local_clean::
	-chmod +w *.depl
	-rm *.depl *.sal tmp_file*.tmp &> /dev/null
