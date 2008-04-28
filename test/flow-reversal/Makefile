#__INSERT_LICENSE__

.PHONY: force all run $(tests) 

default: run

tests: test_right test_left

CASE_NAME := strip2d
NP := 1
MACHI = -machinefile machi.dat

LOC_PROG = ns
PROG = $(NS)

BOPT := g_c++
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test
SRCS = 

mesh:
	make-depl strip2d.epl strip2d.depl -d mkmesh=1

run: 
	mkstepsdir.pl
	make-depl strip2d.epl strip2d.depl -d mkmesh=0
	-$(MPIEXEC) $(NL) -np $(NP) $(PROG) -case $(CASE_NAME).depl 

local_clean::
	-putex
	-rm -rf mat.output system.dat save.state data.m \"\"

link:
	ln -sf $(PROG) $(LOC_PROG)

prog:
	$(MAKE) -C $(PETSCFEM_DIR) ns

#s
