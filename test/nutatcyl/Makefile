#__INSERT_LICENSE__
# $Id: Makefile,v 1.2 2002/08/17 22:39:00 mstorti Exp $ 

.PHONY: force all run nutcyl tests

PETSCFEM_DIR := ../..

tests: nutcyl

outfile = system.dat.tmp

NAME := closcyl
NP := 1

export NAME NP

LOC_PROG = ns_O
PROG = $(PETSCFEM_DIR)/applications/ns/$(LOC_PROG).bin
CASE_NAME = closcyl

patch:
	$(MAKE) -C $(HOME)/SOFT/matlab patch
	cp $(HOME)/SOFT/matlab/README.shar matlab.patch

# PETSCFEM_DIR = $(HOME)/PETSC/petscfem-beta-2.94.pl8
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test
SRCS = 

data: $(NAME).depl force
$(NAME).depl: force

run: data
	-rm -f $(NAME).some_rslt.tmp
	if [ -e $(NAME).force.tmp ] ; then rm $(NAME).force.tmp ; fi
	if [ -e save.state ] ; then cp save.state save.state.tmp ; fi
	-$(MPIRUN) -np $(NP) $(NL) -machinefile machi.dat $(PROG) \
		 -case $(NAME).depl  | tee petscfem.output.tmp

local_clean::
	-putex
	-rm -rf mat.output system.dat save.state data.m \"\"

nutcyl:
	-rm -f cylinder.force.nuta15.tmp
	$(MAKE) nuta=15 mesh=4-10-I run
	cp cylinder.force.tmp cylinder.force.nuta15.tmp
	-rm -f cylinder.force.nuta0.tmp
	$(MAKE) nuta=0 mesh=4-10-I run
	cp cylinder.force.tmp cylinder.force.nuta0.tmp
	octave -q checknutcyl.m > checknutcyl.verif.tmp
#s
