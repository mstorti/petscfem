#__INSERT_LICENSE__
# $Id: Makefile,v 1.2 2001/07/28 20:04:14 mstorti Exp $ 

.PHONY: force tangent stabi data auxi

run: 
	$(MAKE) base_flow flux_jaco 	
	mv system.dat flux_jaco.dat
	$(MAKE) mass_jaco

PETSCFEM_DIR = ../../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test
SRC_DIR = $(PETSCFEM_DIR)/src
PETSCFEM_LIB = $(SRC_DIR)/libpetscfem.a
SRCS = 

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
DATA_FILES = stabi.peri.tmp stabi.nod.tmp stabi.con.tmp 
data:
	octave -q < mkmesh.m 

.PHONY: base_flow flux_jaco mass_jaco
base_flow flux_jaco mass_jaco: 
	$(MAKE) -W stabi.epl EPERLFLAGS="-d case=$@" stabi.depl
	-$(MPIRUN) -np 1 $(NS) -case stabi.depl >$@.out.tmp

stabi: 
	$(MAKE) EPERLFLAGS='-d stabi=1' -W stabi.epl stabi.depl
	-$(MPIRUN) -np 1 $(NS) -case stabi.depl >stabi.out.tmp

local_clean::
	-rm system.dat flux_jaco.dat

auxi_o:
	ln -sf $(PETSCFEM_DIR)/applications/ns/ns.bin ns
	echo set args -case stabi.depl >.gdbinit
	echo file ns >>.gdbinit
