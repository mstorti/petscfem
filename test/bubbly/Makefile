# $Id: Makefile,v 1.9.80.1 2006/02/24 12:26:57 rodrigop Exp $ 

.PHONY: all run  compile lclean laplace adv nso laplace_tests force \
        ctest gprof

tests: test_sqcav

all: run

# PETSCFEM_DIR = $(HOME)/PETSC/petscfem
PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

ADVDIF  = $(PETSCFEM_DIR)/../petscfem/applications/advdif/advdif_g.bin
## ultima version mas o menos estable al momento de hacer la implementacion original  
ADVOLD  = $(PETSCFEM_DIR)/../petscfem.2-94.pl20/applications/advdif/advdif_g.bin

LOC_PROG := advdif_g
#CASE_NAME := bubbly
CASE_NAME := poise2fases
# CASE_NAME := couette
# PROG := $(NS)
PROG := $(ADVOLD)
PROG := $(ADVDIF)


CASE := $(CASE_NAME)

################################################################
NP = 1

MCH =    

NP := $(shell procsel $(PETSCFEM_PROCTABLE))
ifeq ($(np),)
np := 1
endif

#w Runs program
run:  data
	-${MPIRUN} -np $(np) $(MCH) $(PROG) \
			-case $(CASE).depl # > output.out

data: force
	$(MAKE) $(CASE).depl
	octave -q mkpoise2fases.m

# EPERLFLAGS = -d case=$$case
$(CASE).depl: force
#	echo in target \$(CASE).depl : CASE $(CASE)

local_clean::
	-chmod +w *.depl
	-rm *.depl *.sal tmp_file*.tmp &> /dev/null

usens:
	ln -sf mkcouette.m.ns mkcouette.m
	ln -sf couette.epl.ns couette.epl
	ln -sf $(NS) prog

usebb:
	ln -sf mkcouette.m.bb mkcouette.m
	ln -sf couette.epl.bb couette.epl
	ln -sf $(ADVDIF) prog

cycle:
	for ip in 1 2 1 2 1 2 1 2 1 2 ; do export ip ; \
	make run ; done

