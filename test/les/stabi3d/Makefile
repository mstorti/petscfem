#__INSERT_LICENSE__
# $Id: Makefile,v 1.20 2004/07/26 19:54:48 mstorti Exp $ 

.PHONY: force tangent stabi data auxi

run: stabi

CASE := stabi
PETSCFEM_DIR := ../../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

nolocal := yes
include ./Makefile.parproc

stabi: 	
	saverot stabi.some_rslt%d.tmp
	$(MAKE) $(CASE).depl
	$(OCTAVE) < mkmesh.m 
	-$(MPIEXEC) -machinefile ./machi.dat -n $(NP) -nolocal $(PROG) -case stabi.depl

stabi.con-tet.tmp: stabi.con.tmp
	$(MAKE) -C $(PETSCFEM_DIR)/tools/ hexasplit.bin
	$(PETSCFEM_DIR)/tools/hexasplit.bin -xi $< -o $@

dx: viewer_tmp.dx stabi.con-tet.tmp dxdata.tmp

viewer_tmp.dx dxdata.tmp: viewer.dx stabi.state.tmp
	$(EPERL) -P viewer.dx > viewer_tmp.dx

LOC_PROG := ns
CASE_NAME := stabi
PROG := $(PETSCFEM_DIR)/applications/ns/ns_O.bin

local_clean::
	-rm -rf viewer_tmp.dx
