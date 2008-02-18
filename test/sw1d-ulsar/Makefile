#__INSERT_LICENSE__
# $Id: Makefile,v 1.11 2004/09/25 11:32:33 mstorti Exp $ 

default: test_abso

tests: test_abso

CASE=canal1d

# PETSCFEM_DIR := /u/rodrigop/PETSC/petscfem-hydro
PETSCFEM_DIR := ../..

NP := $(shell procsel -n ./proctable)

ADVDIF  = $(PETSCFEM_DIR)/applications/advdif/advdif_O.bin
PROG := $(ADVDIF)
LOC_PROG := advdif

include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

SRCS = 

CLEAN_DIRS = 
DEPEND_DIRS = $(CLEAN_DIRS)

# ADVDIF = $(APPDIR)/advdif/advdif_O.bin

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
# Test for abso shallow water 1d

data:
	make-depl $(CASE).epl $(CASE).depl 

mesh: data
	$(OCTAVE) $(CASE).m

run: test_abso

test_abso:
	mkstepsdir.pl
	make-depl $(CASE).epl $(CASE).depl 
	$(MPI_HOME)/bin/mpiexec  -n 1 $(PROG) -case $(CASE).depl
#	time $(MPI_HOME)/bin/mpirun -1 -machinefile ./machi.dat -n $(NP) $(PROG) -case $(CASE).depl
#	$(OCTAVE)  check_test_abso.m

chiquito:
	$(MAKE) chiquito.depl
	-$(MPIRUN) -np $(NP)  $(ADVDIF) -case chiquito.depl

spiral:
	$(MAKE) spiral.depl
	$(OCTAVE) spiral.m
	-$(MPIRUN) -np $(NP)  $(ADVDIF) -case spiral.depl

corner:
	$(MAKE) corner.depl
	$(OCTAVE) corner.m
	-$(MPIRUN) -np $(NP)  $(ADVDIF) -case corner.depl

igrec:
	$(MAKE) igrec.depl
	$(OCTAVE) igrec.m
	-$(MPIRUN) -np $(NP)  $(ADVDIF) -case igrec.depl

mpdstart:
	let NP1=$(NP)+1 ; \
	$(MPI_HOME)/bin/mpdboot --file=./machi.dat -n $$NP1 --rsh=ssh


local_clean::
	-chmod +w *.depl
	-rm -f STEPS* *.depl *.tmp &> /dev/null

