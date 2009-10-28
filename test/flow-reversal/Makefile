#__INSERT_LICENSE__

.PHONY: force all run $(tests) 

default: run

tests: test_flowrev

CASE_NAME := strip2d
NP := 1
MACHI = -machinefile machi.dat

LOC_PROG = ns
PROG = $(NS)

BOPT := g_c++
PETSCFEM_DIR := ../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test
SRCS = 

test_flowrev:
	make-depl strip2d.epl strip2d.depl -d mkmesh=1 -d gb=1
	-$(MPIEXEC) $(NL) -np $(NP) $(PROG) -case $(CASE_NAME).depl \
		$(call redir,strip2d.outputr.tmp)
	if [ -f strip2d.state.tmp ] ; \
		then cp strip2d.state.tmp strip2d.state-right.tmp ; fi
	make-depl strip2d.epl strip2d.depl -d mkmesh=1 -d gb=-1
	-$(MPIEXEC) $(NL) -np $(NP) $(PROG) -case $(CASE_NAME).depl \
		$(call redir,strip2d.outputl.tmp)
	if [ -f strip2d.state.tmp ] ; \
		then cp strip2d.state.tmp strip2d.state-left.tmp ; fi
	octave -qH check.m > flowrev-verif.tmp

mesh:
	make-depl strip2d.epl strip2d.depl -d mkmesh=1

run: 
	mkstepsdir.pl
	make-depl strip2d.epl strip2d.depl -d mkmesh=0
	-$(MPIEXEC) $(NL) -np $(NP) $(PROG) -case $(CASE_NAME).depl 

local_clean::
	-putex
	-rm -rf STEPS* mat.output system.dat save.state data.m \"\"

link:
	ln -sf $(PROG) $(LOC_PROG)

prog:
	$(MAKE) -C $(PETSCFEM_DIR) ns

#s
