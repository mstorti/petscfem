#__INSERT_LICENSE__
# $Id: Makefile,v 1.26 2003/03/23 20:52:43 mstorti Exp $ 

.PHONY: force all run tests data

all: test_dif_temp_iisd

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include ../Makefile.base.test
SRCS = 

#s Test/Advdif tests
# ADVDIF    = $(PETSCFEM_DIR)/applications/advdif/advdif.bin 

EPERLFLAGS = -d $(CASE)=1
DF = newff.depl
OF = newff.nod.tmp newff.con.tmp 
data:  $(OF) $(DF) #auxi_dat # para que genere .gdbinit y el link

data.m.tmp: data.m.epl force
newff.depl: data.m.tmp force

################################################################
# TESTS FOR THE NEW FLUX FUNCTIONS...

tests: tests_plain tests_iisd

tests_plain: test_dif_temp test_dif_temp_cp			\
	test_adv_dif_temp test_adv_temp_cp			\
	test_pure_adv test_reac_adv_dif_temp_y			\
	test_reac_steady test_reac_dif_temp test_std_ard_x_y	\
	test_full_jacs test_full_jacs_cp test_full_full_jacs	\
	test_full_full_jacs_wf test_full_full_jacs_cp		\
	test_full_full_jacs_adv_cp				\
	test_full_full_jacs_t					\
	test_full_full_jacs_tr test_adv_dif_stdy_bcconv		\
	test_stdy_dif test_check_fields test_check_nodes 

tests_iisd: test_adv_dif_temp_iisd test_adv_temp_cp_iisd	\
	test_dif_temp_iisd					\
	test_full_full_jacs_iisd test_full_jacs_iisd		\
	test_full_jacs_cp_iisd test_pure_adv_iisd		\
	test_reac_adv_dif_temp_y_iisd test_std_ard_x_y_iisd 

# test suite, used when performing particular tests
testsss: test_dif_temp_iisd test_dif_temp_cp_iisd \
	test_adv_dif_temp_iisd test_adv_temp_cp_iisd

test_%: 
	$(MAKE) CASE=$* run

data_%:
	$(MAKE) CASE=$* data

#w Test for per-element properties 
$(OF): newff.m data.m.tmp
	touch fulljacdef.tmp	# because it is needed by ePerl!!
	$(OCTAVE) newff.m >/dev/null

#w makes auxiliary files .gdbinit and link advdif
auxi_dat: advdif .gdbinit

.gdbinit:
	echo file advdif >.gdbinit
	echo set args -case newff.depl >>.gdbinit

advdif: advdif.link

NPP := $(shell procsel $(PETSCFEM_PROCTABLE))
# use NP=2 for all IISD cases
ifeq ($(findstring _iisd,$(CASE)),_iisd)
NP:= 2
endif
ifeq ($(NP),)
NP := 1
endif

run: data
	-rm -f save.state.tmp 
	-$(MPIRUN) -np $(NP) -machinefile machi.dat $(ADVDIF) \
			-case newff.depl > newff.$(CASE).out
	if [ -e save.state.tmp ] ; then \
		$(OCTAVE) verif.m > $(CASE).vrf.tmp ; fi

local_clean::
	-chmod +w newff.data.res
	-rm -f data.m save.state*

LOC_PROG = advdif
PROG = $(ADVDIF)
CASE_NAME = newff

#s
