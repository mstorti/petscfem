#__INSERT_LICENSE__
# $Id: Makefile,v 1.10 2002/11/23 19:26:26 mstorti Exp $ 

.PHONY: force tangent stabi data auxi

run: stabi

CASE := stabi
PETSCFEM_DIR := ../../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

include ./Makefile.parproc

stabi: 	
	saverot stabi.some_rslt%d.tmp
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
