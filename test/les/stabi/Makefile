#__INSERT_LICENSE__
# $Id: Makefile,v 1.1 2001/05/01 21:47:15 mstorti Exp $ 

.PHONY: fstack force les_data

run: stabi

PETSCFEM_DIR = ../../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test
SRC_DIR = $(PETSCFEM_DIR)/src
PETSCFEM_LIB = $(SRC_DIR)/libpetscfem.a
SRCS = 

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
DATA_FILES = stabi.peri.tmp stabi.nod.tmp stabi.con.tmp 
.PHONY: data
data:
	octave -q < mkmesh.m 
	$(MAKE) -W stabi.epl stabi.depl

.PHONY: stabi
stabi: $(DATA_FILES) stabi.depl
	-$(MPIRUN) -np 1 $(NS) -case stabi.depl 

local_clean::
	-rm les.flag 

