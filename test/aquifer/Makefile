# $Id: Makefile,v 1.27 2002/12/02 19:09:45 mstorti Exp $ 

.PHONY:  force data test_plain test_var_eta test_transient rung run runs

SRCS = rain.cpp rain2.cpp fun.cpp fun2.cpp fun3.cpp

#w All tests
tests: test_plain test_var_eta test_transient test_stream \
		test_as test_dl_fun test_var_rain

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
# CASE_NAME := stream
CASE_NAME := aqui
PROG := advdif

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
	octave -qH check.m > plain.verif.tmp
	$(MAKE) casen=var_eta0 run
	octave -qH check.m > var_eta0.verif.tmp
	octave -qH checkc.m > const.verif.tmp

#w case with variable aquifer depth
test_var_eta:
	$(MAKE) casen=var_eta1 run
	octave -qH check.m > var_eta1.verif.tmp

#w Test transient with constant depth
test_transient:
	$(MAKE) casen=transient run
	octave -qH checku.m > transient.verif.tmp

#w Test for the stream elemset
test_stream_chezy_rect: fun.efn fun2.efn fun3.efn
	$(MAKE) casen=chezy_rect runs

#w Test for the stream elemset
test_stream:
	$(MAKE) fric=chezy shape=rect casen=chezy_rect runs
	$(MAKE) fric=manning shape=circular casen=mann_circ runs
	octave -qH checks.m > stream.verif.tmp

#w Test the coupled aquifer+stream test
test_as:
	$(MAKE) casen=aquist runas
	octave -qH checkas.m > aquist.verif.tmp

.PHONY: funs test_dl_fun1 test_dl_fun2 test_dl_fun3
funs: fun.efn fun2.efn fun3.efn

test_dl_fun: test_dl_fun1 test_dl_fun2 test_dl_fun3

#w Test `dl_generic' amplitude function. Test `fun1'
test_dl_fun1: funs 
	$(MAKE) fric=chezy shape=rect casen=fun1 runs
	octave -qH checkdlf.m > output.dl_fun1.verif.tmp

#w Test `dl_generic' amplitude function. Test `fun2'
test_dl_fun2: funs 
	$(MAKE) fric=chezy shape=rect casen=fun2 runs
	octave -qH checkdlf.m > output.dl_fun2.verif.tmp

#w Test `dl_generic' amplitude function.  Test `fun3'
test_dl_fun3: funs 
	$(MAKE) fric=chezy shape=rect casen=fun3 runs
	octave -qH checkdlf.m > output.dl_fun3.verif.tmp

#w case with constant aquifer depth
test_var_rain: rain2.efn
	-rm -f var_rain.some_out.tmp
	$(MAKE) casen=var_rain run
	octave -qH check_var_rain.m > var_rain.verif.tmp

#w Make data
data: force
	$(MAKE) -W $(CASE).epl $(CASE).depl
	octave -qH $(mkm).m

#w Runs program
rung: data force
	-${MPIRUN} -np $(np) $(MCH) $(PROG) \
			-case $(CASE).depl > output.$(casen).out

# EPERLFLAGS = -d case=$$case
$(CASE).depl: force
#	echo in target \$(CASE).depl : CASE $(CASE)

local_clean::
	-chmod +w *.depl
	-rm *.depl *.sal tmp_file*.tmp &> /dev/null
