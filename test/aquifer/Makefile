# $Id: Makefile,v 1.12 2002/02/10 22:33:28 mstorti Exp $ 

.PHONY:  force data test_plain test_var_eta test_transient rung run runs

#w All tests
tests: test_plain test_var_eta test_transient

#w Runs program for the `aquifer' tests
run: CASE=aqui
run: mkm=mkmesh
run: rung

#w Runs program for the `aquifer_stream' tests
runas: CASE=aquist
runas: mkm=mkmeshaqs
runas: rung

#w Runs program for the `stream' tests
runs: CASE=stream
runs: mkm=mkmeshs
runs: rung

LOC_PROG := advdif
CASE_NAME := stream
PROG := $(ADVDIF)

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

################################################################
NP = 1

PROG := $(ADVDIF)
MCH = -machinefile ./machi.dat

NP := $(shell procsel $(PETSCFEM_PROCTABLE))
ifeq ($(np),)
np := 1
endif

#w case with constant aquifer depth
test_plain:
	$(MAKE) casen=plain run
	octave -q check.m > plain.verif.tmp
	$(MAKE) casen=var_eta0 run
	octave -q check.m > var_eta0.verif.tmp
	octave -q checkc.m > const.verif.tmp

#w case with variable aquifer depth
test_var_eta:
	$(MAKE) casen=var_eta1 run
	octave -q check.m > var_eta1.verif.tmp

#w Test transient with constant depth
test_transient:
	$(MAKE) casen=transient run

#w Test for the stream elemset
test_stream_chezy_rect: fun.efn fun2.efn fun3.efn
	$(MAKE) casen=chezy_rect runs

#w Test the coupled aquifer+stream test
test_as:
	$(MAKE) runas

#w Make data
data: force
	$(MAKE) -W $(CASE).epl $(CASE).depl
	octave -q $(mkm).m

#w Runs program
rung: data force
	-${MPIRUN} -np $(np) $(MCH) $(PROG) \
			-case $(CASE).depl # > output.$(casen).out

# EPERLFLAGS = -d case=$$case
$(CASE).depl: force
#	echo in target \$(CASE).depl : CASE $(CASE)

local_clean::
	-chmod +w *.depl
	-rm *.depl *.sal tmp_file*.tmp &> /dev/null

