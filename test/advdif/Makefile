#__INSERT_LICENSE__

.PHONY: plano data run sine 

all: run

tests: sine

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include ../Makefile.base.test
SRCS = 

#s Test/Advdif tests
#ADVDIF    = $(PETSCFEM_DIR)/applications/advdif/advdif.bin 
LOC_PROG := advdif
CASE_NAME := advec
PROG := $(ADVDIF)

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
	$(EPERL) -P $(case).epl > $(case).depl
	-chmod -w $(eperl_out_$(case))

$(octave_out_$(case)): $(octave_in_$(case))
ifneq ($(octave_out_$(case)),)
	-rm -f $(octave_out_$(case))
endif
	$(OCTAVE) $(octave_script_$(case)) >/dev/null
	-chmod -w $(octave_out_$(case))

sine: sine_fine_mesh sine_crank_nic

sine_fine_mesh:
	$(EPERL) -d nx=20 -d ny=10 -d N_step_period=32 sine.data.epl >sine.data
	$(MAKE) case=sine run
	$(OCTAVE) sine_fine_test.m > sine_fine_mesh.out

sine_crank_nic: 
	$(MAKE) N_step_period=16 sine_crank_nic_npp
	$(MAKE) N_step_period=32 sine_crank_nic_npp
	$(MAKE) N_step_period=128 sine_crank_nic_npp
	$(OCTAVE) sine_crank_nic.m >  sine_crank_nic.out

sine_crank_nic_npp:
	$(EPERL) -d nx=10 -d ny=5 -d N_step_period=$(N_step_period) -d tol_error=0 \
			sine.data.epl >sine.data
	-rm -f save.state
	$(MAKE) case=sine run
	mv -fv save.state save.state.$(N_step_period)

# plano.depl: plano.epl plano.m
# burgers.depl: burgers.epl plano.m
# sine.depl: sine.epl sine.m

plano.con plano.nod: plano.m
	$(OCTAVE) plano.m >/dev/null

advec: advec.depl
	-$(MPIEXEC) -n 1 -machinefile machi.dat $(ADVDIF) \
			-case advec.depl $(call in_backg,$@)

# .gdbinit:
# 	echo file advdif >$@
# 	echo set args -case per_elem_prop.depl >> $@

# DF = per_elem_prop.depl
# OF = sine.nod.tmp
# pdata: $(DF) $(OF)

sine.data: sine.data.epl
	$(EPERL) -P -d per_elem_prop=1 $< > $@

# ################################################################
# # TESTS FOR THE NEW FLUX FUNCTIONS...
# # 1./ DRS: Diffusion reaction with source term. All
# #         combinations of possibilities in the representation. 
# # 2./ AD: Advection-diffusion. All
# #         combinations of possibilities in the representation. 
# # 3./ The same, with mesh rotated 45 deg. (non-isotropic diffusion). 
# # 4./ The same with full jacobians obtained rotating
# #        the components. 


# #w Test for per-element properties 
# $(DF): sine_new.epl sine.data
# 	if [ -e $(DF) ] ; then chmod +w $(DF) ; fi
# 	$(EPERL) -P -d perl_elem_prop=1 sine_new.epl > $(DF)
# 	chmod -w $(DF)

# $(OF): sine_new.m sine.data
# 	$(OCTAVE) sine_new.m

# per_elem_prop: .gdbinit advdif.link $(DF) $(OF)
# 	-$(MPIEXEC) -n 1 -machinefile machi.dat $(ADVDIF) \
# 			-case $(DF) > per_elem_prop.out.tmp

# lclean:
# 	-chmod +w sine.depl
# 	-rm $(octave_files) sine.depl mat.output save.state

#include makefile.d

#NP := $(shell procsel $(PETSCFEM_PROCTABLE))
run: data
	-rm -f dummy $(output_files_$(case))
	-$(MPIEXEC) -n 1 $(ADVDIF) -case $(case).depl </dev/null > $(case).out

# Tests to be done:
# * Convergence to analityc solution
# * Quadratic convergence of CN
# * Invariance with respect of number of processors. 
# * Invariance with respect of weak_form

# tests: testanal

# testanal:
# # For Dt
# 	-rm $(case).some sine_anal_test.sal sine.depl
# 	$(MAKE) data
# 	-$(MPIEXEC) -n 1 -machinefile machi.dat $(ADVDIF) \
# 		 -case $(case).depl >/dev/null
# 	cp save.state save.state.anal.tmp
# 	cp sine.some sine.some.anal.tmp
# 	$(OCTAVE) < sine_anal_test.m > sine_anal_test.out

local_clean::
	-chmod +w sine.data.res
	-rm -f $(eperl_out_$(case)) $(octave_out_$(case)) $(output_files_$(case)) \
		sine.some sine.data sine.data.res save.state*

#w Make a link of the binary to this dir (for use with gdb)
link:
	ln -sf $(PROG) $(LOC_PROG)
#s
