#__INSERT_LICENSE__
# $Id: Makefile,v 1.21 2005/02/17 11:50:44 mstorti Exp $ 

tests = nutcyl test_rigid test_strip \
	test_strip_prism test_strip_prism_inv test_strip2d test_strip2d_inv

.PHONY: force all run $(tests) 

PETSCFEM_DIR := ../..

tests: $(tests)

outfile = system.dat.tmp

ifeq ($(CASE_NAME),)
CASE_NAME := closcyl
endif
NP := 1

ifeq ($(NP),1)
MACHI = -machinefile machi.dat
endif

export CASE_NAME NP use_en use_prisms use_en

LOC_PROG = ns
PROG = $(NS)

# PETSCFEM_DIR = $(HOME)/PETSC/petscfem-beta-2.94.pl8
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test
SRCS = 

data: $(CASE_NAME).depl force
$(CASE_NAME).depl: force

run: data
	-rm -f $(CASE_NAME).some_rslt.tmp
	if [ -e $(CASE_NAME).force.tmp ] ; then rm $(CASE_NAME).force.tmp ; fi
	if [ -e save.state ] ; then cp save.state save.state.tmp ; fi
	-$(MPIRUN) -np $(NP) $(NL) $(MACHI) $(PROG) \
		 -case $(CASE_NAME).depl  > petscfem.output.tmp

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
	$(MAKE) CASE_NAME=strip run
	$(OCTAVE) checkstrip.m > strip.verif.tmp

test_strip_prism: 
	-rm -f stripp.force.tmp stripp.verif.tmp
	$(MAKE) CASE_NAME=strip use_prisms=1 use_en=1 run
	$(OCTAVE) checkstrip.m > stripp.verif.tmp

test_strip_prism_inv: 
	-rm -f strippi.force.tmp strippi.verif.tmp
	$(MAKE) CASE_NAME=strip use_prisms=1 run
	$(OCTAVE) checkstrip.m > strippi.verif.tmp

test_strip2d: 
	-rm -f strip2d.force.tmp strip2d.verif.tmp
	$(MAKE) CASE_NAME=strip2d use_en=1 run
	$(OCTAVE) checkstrip2d.m > strip2d.verif.tmp

test_strip2d_inv: 
	-rm -f strip2di.force.tmp strip2di.verif.tmp
	$(MAKE) CASE_NAME=strip2d use_en=0 run
	$(OCTAVE) checkstrip2d.m > strip2di.verif.tmp

figure:
	$(MAKE) BOPT=O_c++ Omega=3000 Omega_nut=500 nu=0.0572 nuta=20 mesh=10-30-I run

link:
	ln -f $(PROG) $(LOC_PROG)

#s
