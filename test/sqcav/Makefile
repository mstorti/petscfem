# $Id: Makefile,v 1.28 2002/09/07 20:41:23 mstorti Exp $ 

.PHONY: all run  compile lclean laplace adv nso laplace_tests force \
        ctest gprof

runn: 
	$(MAKE) case=weak_form_1 runs

tests: test_re1000 test_iisd test_iisd2 test_meas_perf test_g_body \
		test_fs

PETSCFEM_DIR = $(HOME)/PETSC/petscfem
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

LOC_PROG := ns
#CASE_NAME := qharm
CASE_NAME := sqcav
PROG := $(NS)

################################################################
NP = 1

PROG := $(NS)
MCH = -machinefile ./machi.dat

NP := $(shell procsel $(PETSCFEM_PROCTABLE))
ifeq ($(np),)
np := 1
endif

#w Runs program
run:  data
	-${MPIRUN} -np $(np) $(MCH) $(PROG) -case $(CASE).depl \
	> output.CASE_$(CASE).np_$(np).case_$(case)$(meas_perf).out.tmp

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

test_g_body:
	$(MAKE) case=g_body runs
	octave -q < check.g_body.m > check.g_body.verif.tmp

# tests the square cavity in a 20x20 mesh at Re=1000
# with respect with a reference solution obtained also with Petscfem and
# verified with Ghia results
#w IISD for square cavity test in 20x20 mesh on one processor
test_re1000:
	$(MAKE) case=weak_form_0 runs
	$(MAKE) case=weak_form_1 runs
	octave -q < check_re1000.m > check.re1000.verif.tmp

#w IISD for square cavity test in 20x20 mesh on two processors
test_iisd2: 
	$(MAKE) case=lu np=1 runs
	$(MAKE) case=iisd_sbp1 runs
	$(MAKE) case=iisd_sbp4 runs
	$(MAKE) case=iisd_sbp16 runs
	octave -q < check_iisd2.m > check.iisd.verif.np2.tmp

test_meas_perf:
	$(MAKE) case=lu np=1 meas_perf=1 runs

#w Tests IISD solver on the square cavity problem
test_iisd:
	$(MAKE) case=lu np=1 runs
	$(MAKE) case=iisd_sbp1 runs
	$(MAKE) case=iisd_sbp2 runs
	$(MAKE) case=iisd_sbp4 runs
	$(MAKE) case=iisd_sbp8 runs
	$(MAKE) case=iisd_sbp16 runs
	octave -q < check_iisd.m > check.iisd.verif.np1.tmp

#w Tests Fractional step solver
test_fs:
	$(MAKE) case=fractional_step runs
	octave -q < check_re1000_fs.m > check.re_1000_fs.verif.tmp

#w Makes local link for running with debugger
link:
	ln -sf $(PROG) ns

data: force
	$(MAKE) $(CASE).depl
	octave -q mk$(CASE).m

$(CASE).depl: force

local_clean::
	-chmod +w *.depl
	-rm *.depl *.sal tmp_file*.tmp &> /dev/null
