# $Id: Makefile,v 1.10 2001/01/11 15:12:11 mstorti Exp $ 

.PHONY: plano data run

all: run

tests: sine

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include ../Makefile.base.test
SRCS = 

# with this it doesn't work
ADVDIF    = $(PETSCFEM_DIR)/applications/advdif/advdif.bin 

#case := sine
#case := burgers

#octave_files = sine.con.tmp sine.nod.tmp sine.fixa.tmp sine.bcconv.tmp
#octave_script = $(case).m
#data_file = $(case).depl $(case).data.res

data: data_$(case)

.PHONY: data_$(case)

output_files_$(case) := save.state outvector*.sal

################################################################
eperl_in_sine := sine.epl sine.data
eperl_out_sine := sine.depl sine.data.res
eperl_script_sine := sine.epl

octave_in_sine := sine.data sine.data.res
octave_out_sine := sine.con.tmp sine.nod.tmp sine.fixa.tmp sine.bcconv.tmp 
octave_script_sine := sine.m

output_files_sine := $(output_files_sine) sine.some 
################################################################
eperl_in_burgers := burgers.epl
eperl_out_burgers := burgers.depl
eperl_script_burgers := burgers.epl

octave_in_burgers :=
octave_out_burgers :=
octave_script_burgers :=
################################################################

data_$(case): $(eperl_out_$(case)) $(octave_out_$(case))

#data_burgers: burgers.depl plano.con plano.nod plano.m
#data: plano.depl plano.con plano.nod

$(eperl_out_$(case)): $(eperl_in_$(case))
ifneq ($(eperl_out_$(case)),)
	-rm -f $(eperl_out_$(case))
endif
	eperl -P $(case).epl > $(case).depl
	-chmod -w $(eperl_out_$(case))

$(octave_out_$(case)): $(octave_in_$(case))
ifneq ($(eperl_out_$(case)),)
	-rm -f $(octave_out_$(case))
endif
	octave -q $(octave_script_$(case))
	-chmod -w $(octave_out_$(case))

sine: sine_fine_mesh sine_crank_nic

sine_fine_mesh:
	eperl -d nx=20 -d ny=10 -d N_step_period=32 sine.data.epl >sine.data
	$(MAKE) case=sine run
	octave -q sine_fine_test.m > sine_fine_mesh.out

sine_crank_nic: 
	$(MAKE) N_step_period=16 sine_crank_nic_npp
	$(MAKE) N_step_period=32 sine_crank_nic_npp
	$(MAKE) N_step_period=128 sine_crank_nic_npp
	octave -q sine_crank_nic.m >  sine_crank_nic.out

sine_crank_nic_npp:
	eperl -d nx=10 -d ny=5 -d N_step_period=$(N_step_period) -d tol_error=0 \
			sine.data.epl >sine.data
	$(MAKE) case=sine run
	mv -fv save.state save.state.$(N_step_period)

# plano.depl: plano.epl plano.m
# burgers.depl: burgers.epl plano.m
# sine.depl: sine.epl sine.m

plano.con plano.nod: plano.m
	octave -q plano.m

# lclean:
# 	-chmod +w sine.depl
# 	-rm $(octave_files) sine.depl mat.output save.state

#include makefile.d

NP := $(shell procsel $(PETSCFEM_PROCTABLE))
run: data
	-rm -f dummy $(output_files_$(case))
	-$(MPIRUN) -np 1 -machinefile machi.dat $(ADVDIF) -case $(case).depl > $(case).out

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
	octave -q < sine_anal_test.m > sine_anal_test.out

local_clean::
	-chmod +w sine.data.res
	-rm -f $(eperl_out_$(case)) $(octave_out_$(case)) $(output_files_$(case)) \
		sine.some sine.data sine.data.res save.state*
