# $Id: Makefile,v 1.5 2001/01/04 20:08:31 mstorti Exp $ 

.PHONY: plano data 

all: run

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include ../Makefile.base.test
SRCS = 

# with this it doesn't work
ADVDIF    = $(PETSCFEM_DIR)/applications/advdif/advdif.bin 

case = sine
octave_files = sine.con.tmp sine.nod.tmp sine.fixa.tmp sine.bcconv.tmp
octave_script = $(case).m
data_file = $(case).depl $(case).data.res

data: $(data_file) $(octave_files)
#data: burgers.depl plano.con plano.nod plano.m
#data: plano.depl plano.con plano.nod

$(data_file): $(case).epl $(case).data
	-if [ -e $@ ] ; then chmod +w $(data_file) ; rm $(data_file) ; fi
	eperl -P $< > $@
	-chmod -w $(data_file)

$(octave_files): $(case).m  $(case).data
	octave -q sine.m

plano.depl: plano.epl plano.m
burgers.depl: burgers.epl plano.m
sine.depl: sine.epl sine.m

plano.con plano.nod: plano.m
	octave -q plano.m

# lclean:
# 	-chmod +w sine.depl
# 	-rm $(octave_files) sine.depl mat.output save.state

#include makefile.d

NP := $(shell procsel $(PETSCFEM_PROCTABLE))
run: data
	-rm $(case).some
	-$(MPIRUN) -np 1 -machinefile machi.dat $(ADVDIF) -case $(case).depl

# Tests to be done:
# * Convergence to analityc solution
# * Quadratic vonvergence of CN
# * Invariance with respect of number of processors. 
# * Invariane with respect of weak_form

tests: testanal

testanal:
# For Dt
	-rm $(case).some sine_anal_test.sal sine.depl
	$(MAKE) data
	-$(MPIRUN) -np 1 -machinefile machi.dat $(ADVDIF) \
		 -case $(case).depl >/dev/null
	cp save.state save.state.anal.tmp
	cp sine.some sine.some.anal.tmp
	octave -q < sine_anal_test.m > sine_anal_test.sal

local_clean::
	-chmod +w sine.data.res
	-rm sine.data.res sine.ini 

