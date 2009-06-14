#__INSERT_LICENSE__

.PHONY: force all run $(tests) 

default: tests

tests: test_strip3d

CASE := strip3d
NP := 1
MACHI = -machinefile machi.dat

LOC_PROG = ns
PROG = $(NS)

PETSCFEM_DIR := ../../..

include $(PETSCFEM_DIR)/Makefile.base
MPI_HOME := /usr
include $(PETSCFEM_DIR)/test/Makefile.base.test
SRCS = 

mesh:
	make-depl strip3d.epl strip3d.depl -d mkmesh=1

test_strip3d: 
	make-depl strip3d.epl strip3d.depl -d mkmesh=1
	$(MPIEXEC) $(NL) -np $(NP) $(PROG) -case $(CASE).depl 
# $(call redir,$(CASE).output.tmp)
#	$(OCTAVE) verif.m > embgath-verif.tmp

local_clean::
	-putex
	-rm -rf mat.output system.dat save.state data.m \"\"

link:
	ln -sf $(PROG) $(LOC_PROG)

prog:
	$(MAKE) -C $(PETSCFEM_DIR) ns

#s
