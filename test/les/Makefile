#__INSERT_LICENSE__
# $Id: Makefile,v 1.7 2002/12/02 19:09:46 mstorti Exp $ 

.PHONY: fstack force les_data

run: les

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include ../Makefile.base.test
SRC_DIR = $(PETSCFEM_DIR)/src
PETSCFEM_LIB = $(SRC_DIR)/libpetscfem.a
SRCS = 

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
LES_DATA_FILES = les.nod.tmp  les.con.tmp les.fixa.tmp \
		les.peri.tmp les.wall.tmp les.ini.tmp

MCH = -machinefile ./machi.dat
NP := $(shell procsel $(PETSCFEM_PROCTABLE))
NSO = $(PETSCFEM_DIR)/../petscfem-2.16g/applications/ns/ns.bin
les: les_data
	-$(MPIRUN) -np 1 $(MCH)  $(NS) -log_info -case les.depl | tee output.out
#	-$(MPIRUN) -np $(NP) $(MCH) $(NS) -case les.depl   | tee output.out

prof: pfemprof.txt

pfemprof.txt: gmon.out
#	gprof $(NS) -l gmon.out >$@
	gprof $(NS) gmon.out >$@

les_data: les.depl les.nod.tmp

$(LES_DATA_FILES): les.m
	octave -qH <les.m 
#	touch les.flag

local_clean::
	-rm les.flag 

