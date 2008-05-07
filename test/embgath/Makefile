#__INSERT_LICENSE__

.PHONY: force all run $(tests) 

default: run

tests: test_strip_perelem

CASE_NAME := strip2d
NP := 2
MACHI = -machinefile machi.dat

LOC_PROG = ns
PROG = $(NS)

PETSCFEM_DIR := ../..

BOPT := g_c++
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test
SRCS = 

mesh:
	make-depl strip2d.epl strip2d.depl -d mkmesh=1

test_strip_perelem: 
	make-depl strip2d.epl strip2d.depl -d mkmesh=1
	-$(MPIEXEC) $(NL) -np $(NP) $(PROG) -case $(CASE_NAME).depl \
		$(call redir,strip2d.output.tmp)
	$(OCTAVE) verif.m > embgath-verif.tmp

local_clean::
	-putex
	-rm -rf mat.output system.dat save.state data.m \"\"

link:
	ln -sf $(PROG) $(LOC_PROG)

prog:
	$(MAKE) -C $(PETSCFEM_DIR) ns

#s
