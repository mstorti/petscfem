#__INSERT_LICENSE__
# $Id: Makefile,v 1.13 2003/01/08 15:54:26 mstorti Exp $ 

.PHONY: force all run nutcyl tests

PETSCFEM_DIR := ../..

tests: nutcyl test_rigid test_strip

outfile = system.dat.tmp

NAME := closcyl
NP := 1

ifeq ($(NP),1)
MACHI = -machinefile machi.dat
endif

export NAME NP

LOC_PROG = ns$(osfx)
PROG = $(PETSCFEM_DIR)/applications/ns/$(LOC_PROG).bin
CASE_NAME = strip

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
	-$(MPIRUN) -np $(NP) $(NL) $(MACHI) $(PROG) \
		 -case $(NAME).depl > petscfem.output.tmp

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
	$(OCTAVE) checknutcyl.m > checknutcyl.verif.tmp

test_rigid:
	-rm -f cylinder.force.rigid.tmp rigid.verif.tmp
	$(MAKE) rigid_body_test=1 mesh=4-10-I run
	cp cylinder.force.tmp cylinder.force.rigid.tmp
	$(OCTAVE) checkrigid.m > rigid.verif.tmp

test_strip: 
	-rm -f strip.force.tmp strip.verif.tmp
	$(MAKE) NAME=strip run
	$(OCTAVE) checkstrip.m > strip.verif.tmp

figure:
	$(MAKE) BOPT=O_c++ Omega=3000 Omega_nut=500 nu=0.0572 nuta=20 mesh=10-30-I run

#s
