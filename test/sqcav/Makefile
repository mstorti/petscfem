# $Id: Makefile,v 1.42 2003/02/07 23:18:32 mstorti Exp $ 

.PHONY: all run  compile lclean laplace adv nso laplace_tests force \
        ctest gprof dxhook

dxhook:
	$(MAKE) case=weak_form_0 runs

runn: 
	$(MAKE) case=weak_form_1 runs

tests: test_re1000 test_iisd test_iisd2 test_meas_perf test_g_body test_fs \
	test_zero_w_proc test_zero_w_proc2


# PETSCFEM_DIR = $(HOME)/PETSC/petscfem
PETSCFEM_DIR = ../..
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
	 # > output.CASE_$(CASE).np_$(np).case_$(case)$(meas_perf).out.tmp

export CASE
runs: CASE := sqcav
runs:
	$(MAKE) run

runs2: CASE := sqcav2
runs2:
	$(MAKE) run

runq: CASE := qharm
runq: force
#	$(MAKE) case=iisd_sbp1 run
	for n in 1 2 4 6 8 16 32 64 128 ; \
		do time $(MAKE) case=iisd_sbp$$n run ; \
		mv output.out output.iisd$$n.out ; \
	done

runq1: CASE := qharm
runq1: 
	$(MAKE) case=iisd_sbp1 run

runs: CASE := sqcav
test_big:
	time $(MAKE) case=big runs

test_g_body:
	$(MAKE) case=g_body runs
	$(OCTAVE) < check.g_body.m > check.g_body.verif.tmp

# tests the square cavity in a 20x20 mesh at Re=1000
# with respect with a reference solution obtained also with Petscfem and
# verified with Ghia results
#w IISD for square cavity test in 20x20 mesh on one processor
test_re1000:
	$(MAKE) case=weak_form_0 runs
	$(MAKE) case=weak_form_1 runs
	$(OCTAVE) < check_re1000.m > check.re1000.verif.tmp

#w IISD for square cavity test in 20x20 mesh on two processors
test_iisd2: 
	$(MAKE) case=lu np=1 runs
	$(MAKE) case=iisd_sbp1 runs
	$(MAKE) case=iisd_sbp4 runs
	$(MAKE) case=iisd_sbp16 runs
	$(OCTAVE) < check_iisd2.m > check.iisd.verif.np2.tmp

# REDIR := >
REDIR := | tee

#w Tests processors with zero eqs./elements (zero weight)
test_zero_w_proc: CASE := sqcav
test_zero_w_proc:
# 	Runs with one and two process and compares
#
#	One process:
	$(MAKE) case=zwproc_ref np=1 data
	echo $(HOSTNAME) > machi.dat
	echo 1. > weights.dat
	-${MPIRUN} -np 1 -machinefile ./machi.dat $(PROG) -case $(CASE).depl \
	$(REDIR) output.zwproc_ref.out.tmp
#	Two process:
	$(MAKE) case=zwproc np=2 data
	echo $(HOSTNAME) > machi.dat
	echo 0. > weights.dat
	echo 1. >> weights.dat
	-${MPIRUN} -np 2 -machinefile ./machi.dat $(PROG) -case $(CASE).depl \
	$(REDIR) output.zwproc.out.tmp
#
	$(OCTAVE) check_zwproc.m > check.zwproc.tmp

#w Tests processors with zero eqs./elements (zero weight)
test_zero_w_proc2: CASE := sqcav
test_zero_w_proc2:
# 	Runs with one and two process and compares
#
#	One process:
	$(MAKE) case=zwproc2_ref np=1 data
	echo $(HOSTNAME) > machi.dat
	echo 1. >> weights.dat
	-${MPIRUN} -np 1 -machinefile ./machi.dat $(PROG) -case $(CASE).depl \
	$(REDIR) output.zwproc2_ref.out.tmp
#
#	Two process:
	$(MAKE) case=zwproc2 np=2 data
	echo $(HOSTNAME) > machi.dat
	echo 0. > weights.dat
	echo 1. >> weights.dat
	-${MPIRUN} -np 2 -machinefile ./machi.dat $(PROG) -case $(CASE).depl \
	$(REDIR) output.zwproc2.out.tmp
#
	$(OCTAVE) check_zwproc2.m > check.zwproc2.tmp

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
	$(OCTAVE) < check_iisd.m > check.iisd.verif.np1.tmp

#w Tests square cavity at Re=400 with fractional step solver. 
test_fs:
	make case=fractional_step_re400 runs2
	$(OCTAVE) < check_fs.m > check.fs.verif.tmp

#w Makes local link for running with debugger
link:
	ln -sf $(PROG) ns

data: force
	$(MAKE) $(CASE).depl
	$(OCTAVE) mk$(CASE).m

# EPERLFLAGS = -d case=$$case
$(CASE).depl: force
#	echo in target \$(CASE).depl : CASE $(CASE)

local_clean::
	-chmod +w *.depl
	-rm *.depl *.sal tmp_file*.tmp &> /dev/null
