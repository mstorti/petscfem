#__INSERT_LICENSE__
# $Id: Makefile,v 1.7 2002/11/22 12:07:34 mstorti Exp $ 

.PHONY: force tangent stabi data auxi

run: stabi

CASE := stabi
PETSCFEM_DIR = /u/mstorti/PETSC/petscfem.mario
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

NP := $(shell procsel -n $(PETSCFEM_PROCTABLE))

stabi: 	
	$(MAKE) $(CASE).depl
	octave -q < mkmesh.m 
	-$(MPIRUN) -machinefile ./machi.dat -np $(NP) -nolocal $(PROG) -case stabi.depl

stabi.con-tet.tmp: stabi.con.tmp
	$(PETSCFEM_DIR)/tools/hexasplit.bin -xi $< -o $@

dx: viewer_tmp.dx

viewer_tmp.dx: stabi.con-tet.tmp stabi.state.tmp viewer.dx
	eperl -P viewer.dx > $@

LOC_PROG := ns
CASE_NAME := stabi
PROG := $(PETSCFEM_DIR)/applications/ns/ns_O.bin

local_clean::
	-rm -rf viewer_tmp.dx
