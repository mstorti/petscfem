#__INSERT_LICENSE__
# $Id: Makefile,v 1.3 2001/05/02 01:42:16 mstorti Exp $ 

.PHONY: force tangent stabi data

run: tangent stabi

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

tangent: 
	$(MAKE) -W stabi.epl stabi.depl
	-$(MPIRUN) -np 1 $(NS) -case stabi.depl 

stabi: 
	$(MAKE) EPERLFLAGS='-d stabi=1' -W stabi.epl stabi.depl
	-$(MPIRUN) -np 1 $(NS) -case stabi.depl 

local_clean::
	-rm system.dat

