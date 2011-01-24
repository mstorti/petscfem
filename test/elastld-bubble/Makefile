#$Id: Makefile,v 1.1 2006/03/11 23:45:43 mstorti Exp $

.PHONY: force all run run_struct

default: run

tests: test_bubble

## activar para compilar la opcion sin debager
# BOPT := g_c++
PETSCFEM_DIR := ../..
SRCS =
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test
MPI_HOME := /usr

CASE := bubble
PROG := $(NS)
LOC_PROG := ns

ifeq ($(HOSTNAME),aquiles)
NP := $(shell procsel -n ./proctable)
MCH := -1 -machinefile ./machi.dat
else
NP := 1
MCH := 
endif

test_bubble:
	$(MAKE) mesh
	$(MAKE) run $(call redir,$@.output.tmp)
	$(OCTAVE) verif.m > $@.check.tmp

mesh:
	make-depl $(CASE).epl $(CASE).depl -dmkmesh=1

split:
	$(PETSCFEM_DIR)/tools/hexasplit.bin \
		-i bubble.con.tmp -o bubble.con-tetra.tmp

run: 
	make-depl $(CASE).epl $(CASE).depl -dmkmesh=0
	$(MPIEXEC) -a $(CASE) $(MCH) -np $(NP) $(PROG) \
		-code struct -case $(CASE).depl

zip:
	shopt -s nullglob ;				\
	find STEPS* -name '$(CASE)*.state*.tmp' |	\
		xargs -r gzip -f 


shell_init:

shell_pre:

shell_post: zip

shell_close:

link:
	$(MAKE) -C $(PETSCFEM_DIR) ns
	ln -sf $(NS) ns

kill:
	-mpdkilljob -a $(CASE) &> /dev/null

ltag:
	makeltag save-log.txt notes.txt

local_clean::
	shopt -s nullglob ; rm -f *.fifo ns advdif *.log *.vtk
	rrm -rf STEPS*

video:
	mkavi.sh bubble

prog:
	$(MAKE) -C $(PETSCFEM_DIR) ns
