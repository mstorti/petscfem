#__INSERT_LICENSE__
# $Id: Makefile,v 1.16 2003/03/19 21:28:07 mstorti Exp $ 

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
#NP := $(shell procsel ./proctable)
PROG := $(NS)

les: les_data
#	-$(MPIRUN) -np $(NP) $(MCH)  $(PROG) -case les.depl | tee output.out
	-$(MPIRUN) -np 2 $(MCH)  $(PROG) -case les.depl | tee output.out

prof: pfemprof.txt

pfemprof.txt: gmon.out
#	gprof $(NS) -l gmon.out >$@
	gprof $(NS) gmon.out >$@

les_data: les.depl les.nod.tmp

$(LES_DATA_FILES): les.m
	$(OCTAVE) <les.m 
#	touch les.flag

local_clean::
	-rm les.flag 

