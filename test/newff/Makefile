# $Id: Makefile,v 1.3 2001/03/28 16:34:28 mstorti Exp $ 

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
data: $(DF) $(OF) #auxi_dat # para que genere .gdbinit y el link

data.m: data.m.epl force
newff.depl: data.m force

################################################################
# TESTS FOR THE NEW FLUX FUNCTIONS...
# 1./ DRS: Diffusion reaction with source term. All
#         combinations of possibilities in the representation. 
# 2./ AD: Advection-diffusion. All
#         combinations of possibilities in the representation. 
# 3./ The same, with mesh rotated 45 deg. (non-isotropic diffusion). 
# 4./ The same with full jacobians obtained rotating
#        the components. 

tests: test_dif_temp test_adv_dif_temp test_reac_adv_dif_temp_y \
	test_reac_steady test_reac_dif_temp test_std_ard_x_y

test_%: 
	$(MAKE) CASE=$* run

#w Test for per-element properties 
$(OF): newff.m data.m
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
	octave -q verif.m >$(CASE).vrf.tmp

local_clean::
	-chmod +w newff.data.res
	-rm -f data.m save.state*

#s
