#__INSERT_LICENSE__
# $Id: Makefile,v 1.13 2001/04/07 21:44:04 mstorti Exp $ 

.PHONY: force all run tests data

all: run

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include ../Makefile.base.test
SRCS = 

#s Test/Advdif tests
ADVDIF    = $(PETSCFEM_DIR)/applications/advdif/advdif.bin 

EPERLFLAGS = -d $(CASE)=1
DF = newff.depl
OF = newff.nod.tmp newff.con.tmp 
data:  $(OF) $(DF) #auxi_dat # para que genere .gdbinit y el link

data.m: data.m.epl force
newff.depl: data.m force

################################################################
# TESTS FOR THE NEW FLUX FUNCTIONS...

tests: test_dif_temp test_dif_temp_cp \
       test_adv_dif_temp test_pure_adv test_reac_adv_dif_temp_y \
	test_reac_steady test_reac_dif_temp test_std_ard_x_y \
	test_full_jacs test_full_jacs_cp test_full_full_jacs \
	test_full_full_jacs_wf test_full_full_jacs_cp \
	test_full_full_jacs_t \
	test_full_full_jacs_tr test_adv_dif_stdy_bcconv \
	test_stdy_dif

testsss: test_full_full_jacs_wf test_full_full_jacs_cp \
	test_full_full_jacs_t \
	test_full_full_jacs_tr test_adv_dif_stdy_bcconv \
	test_stdy_dif

test_%: 
	$(MAKE) CASE=$* run

#w Test for per-element properties 
$(OF): newff.m data.m
	touch fulljacdef.tmp	# because it is needed by eperl!!
	octave -q newff.m

#w makes auxiliary files .gdbinit and link advdif
auxi_dat: advdif .gdbinit

.gdbinit:
	echo file advdif >.gdbinit
	echo set args -case newff.depl >>.gdbinit

advdif: advdif.link

NP := $(shell procsel $(PETSCFEM_PROCTABLE))
run: data
	-$(MPIRUN) -np 1 -machinefile machi.dat $(ADVDIF) \
			-case newff.depl > newff.out
	octave -q verif.m | tee  $(CASE).vrf.tmp

local_clean::
	-chmod +w newff.data.res
	-rm -f data.m save.state*

#s
