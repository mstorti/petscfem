# $Id: Makefile,v 1.2 2002/01/16 01:19:06 mstorti Exp $ 

.PHONY:  force test_plain test_var_eta

#w All tests
tests: test_plain test_var_eta

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

LOC_PROG := advdif
CASE_NAME := aqui
PROG := $(ADVDIF)

CASE = $(CASE_NAME)
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

#w Runs program
run:  data
	-${MPIRUN} -np $(np) $(MCH) $(PROG) -case $(CASE).depl > output.out

data: force
	$(MAKE) $(CASE).depl
	octave -q mkmesh.m

# EPERLFLAGS = -d case=$$case
$(CASE).depl: force
#	echo in target \$(CASE).depl : CASE $(CASE)

local_clean::
	-chmod +w *.depl
	-rm *.depl *.sal tmp_file*.tmp &> /dev/null
