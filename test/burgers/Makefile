#__INSERT_LICENSE__

.PHONY: force all run tests data

default: fd_jac_3

all: run

test_%:
	$(MAKE) CASE=$* run

#w tests: test__d05 test__d01 test__d001
tests: test__d01 test__d001 fd_jac jaco_tests

#w Burger test with viscosity 0.1
test__d01: test__d01_wf test__d01_nwf 
	CASE=_d01 $(OCTAVE) verif.m > $@.out.tmp

#w Burger test with viscosity 0.5
test__d05: test__d05_wf test__d05_nwf 
	CASE=_d05 $(OCTAVE) verif.m > $@.out.tmp

#w Burger test with viscosity 0.01
test__d001: test__d001_wf test__d001_nwf 
	CASE=_d001 $(OCTAVE) verif.m > $@.out.tmp

#w All tests for the  `compute_fd_adv_jacobian' feature
fd_jac: fd_jac_1 fd_jac_2 fd_jac_3

#w Test `compute_fd_adv_jacobian' feature with
#w verbose printing and threshold
fd_jac_1:
	$(MAKE) CASE=fd_jac_1_ref run
	$(MAKE) CASE=fd_jac_1 run
	octave -qH verif_fd_jac_1.m > fd_jac_1.check.tmp

#w Test `compute_fd_adv_jacobian' feature with
#w printing of jacobians and higher threshold
fd_jac_2:
	$(MAKE) CASE=fd_jac_2 fd_jac=2 run

#w Test `compute_fd_adv_jacobian' feature with
#w random elements
fd_jac_3:
	$(MAKE) CASE=fd_jac_3 fd_jac=3 run

#---<*>---//---<*>---//---<*>---//---<*>---//---<*>---// 

# BOPT := g_c++
PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include ../Makefile.base.test
SRCS = 

#s Test/Advdif tests
# ADVDIF    = $(PETSCFEM_DIR)/applications/advdif/advdif.bin 

OF = burgers.nod.tmp burgers.con.tmp burgers.ini.tmp
$(OF): mmesh.m data.m

data: burgers.depl $(OF) 
burgers.depl: burgers.epl data.m force

################################################################
# TESTS FOR THE NEW FLUX FUNCTIONS...

# for use in background
IN_BACKG := </dev/null >burgers.$(CASE).out.tmp
# for debugging
# IN_BACKG :=

NP := $(shell procsel $(PETSCFEM_PROCTABLE))
run: 
	make-depl burgers.epl burgers.depl 
	-$(MPIEXEC) -n 1 $(ADVDIF) -case burgers.depl $(IN_BACKG)
	if [ -f save.state.tmp ] ;					\
		then mv -f save.state.tmp save.state.$(CASE).tmp ;	\
	fi

jaco_tests: 
	$(MAKE) -C jacotest tests

CLEAN_DIRS := jacotest

PROG := $(ADVDIF)
LOC_PROG := advdif

prog:
	$(MAKE) -C $(PETSCFEM_DIR) BOPT=g_c++ $(LOC_PROG)

#w Makes local link for running with debugger
link:
	$(MAKE) -C $(PETSCFEM_DIR) $(LOC_PROG)
	ln -f $(PROG) $(LOC_PROG)

#s
