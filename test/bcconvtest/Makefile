#__INSERT_LICENSE__
# $Id: Makefile,v 1.1 2002/10/09 13:38:44 mstorti Exp $ 

.PHONY: force all run runo tests data

outfile = system.dat.tmp
nomega = 10
export nomega
#jomega = undef
export jomega
runf:
	$(MAKE) NAME=sloshf PROG=$(NSFREQ) MKMESH=mkmeshf.m NP=1 test_case
	octave -q procf.m 

runfa:
	for jomega in `seq 1 $$nomega` ; do \
	echo jomega: $$jomega ; \
	$(MAKE) NAME=sloshf PROG=$(NSFREQ) MKMESH=mkmeshf.m NP=1 test_case ; \
	grep omega sloshf.depl ; \
	octave -q procf.m ; \
	mv sloshfa.rslt.tmp sloshfa.$$jomega.tmp ; \
	mv sloshf.state.tmp sloshfa.state.$$jomega.tmp ; \
	done

runfo:
	$(MAKE) NAME=freq PROG=$(NSFREQ) MKMESH=mkmshf.m NP=1 test_case

#NSFREQ = ../SLOSHING/rockns/rockns.bin
PROG = $(NS)
#NAME = bcconv_test_3d
#MKMESH = mkmesh.m
NAME = bcconv_test_3d
MKMESH = mkmesh.m

run:
#	$(MAKE) NAME=$(NAME) PROG=$(NS) MKMESH=mkmesh.m NP=1 test_case
#	$(MAKE) NAME=$(NAME) PROG=$(NS) MKMESH=mksphere_axi.m NP=1 test_case
# 	$(MAKE) NAME=$(NAME) PROG=$(NS) MKMESH=mksphere_axi_rot.m NP=1 test_case
	$(MAKE) NAME=$(NAME) PROG=$(NS) MKMESH=$(MKMESH) NP=1 test_case

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test
SRCS = 
LOC_PROG = ns
CASE_NAME = $(NAME)

#PROG = ../SLOSHING/rockns/rockns.bin
#LOC_PROG = rockns
CLEAN_DIRS = 

.PHONY: data force
test_case_data: 
#	-if [ ! "$(cont)" = "y" ] ; then rm -f *.tmp ; fi 
	$(MAKE) data
	octave -q $(MKMESH) > /dev/null

data: $(NAME).depl force
$(NAME).depl: force

test_case: test_case_data
	-rm -f slosh.some_rslt.tmp
	if [ -e save.state ] ; then cp save.state save.state.tmp ; fi
	-$(MPIRUN) -np $(NP) -machinefile machi.dat $(PROG) \
			-case $(NAME).depl # | tee output..tmp

local_clean::
	-putex
	-rm -rf mat.output system.dat save.state data.m 
auxia:
	if [ ! -L anima ] ; \
	 then ln -sf ../SLOSHING/mesa/anima.bin anima ; \
	fi
	echo file anima >.gdbinit
	echo set args sloshf.drv >>.gdbinit


#s
