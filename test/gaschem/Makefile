# $Id: Makefile,v 1.12.58.1 2006/02/24 12:26:57 rodrigop Exp $ 

.PHONY: run default

default: run

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

LOC_PROG := advdif
CASE := pool
PROG := $(ADVDIF)

#---<*>---//---<*>---//---<*>---//---<*>---//---<*>---// 
MCH =    

# PETSCFEM_PROCTABLE := ./proctable

#---<*>---//---<*>---//---<*>---//---<*>---//---<*>---// 
ifeq (1,1)
## Sequential run
NP := $(shell procsel -l $(PETSCFEM_PROCTABLE))
NOLOCAL := 
NPMAX := 1
endif
#---<*>---//---<*>---//---<*>---//---<*>---//---<*>---// 
ifeq (1,0)
## Parallell run `no-local'
NP := $(shell procsel -n $(PETSCFEM_PROCTABLE))
NOLOCAL :=  
endif
#---<*>---//---<*>---//---<*>---//---<*>---//---<*>---// 
## Parallell run `local'
ifeq (1,0)
NP := $(shell procsel -l $(PETSCFEM_PROCTABLE))
NOLOCAL := 
endif
#---<*>---//---<*>---//---<*>---//---<*>---//---<*>---// 

NP := $(shell perl -e "print $(NP)> $(NPMAX) ? $(NPMAX) : $(NP);")

#w Runs program
run:
	-rrm -f pool.state_*.tmp
	make-depl $(CASE).epl $(CASE).depl
	octave -qH mkmesh.m
	-${MPIRUN} $(NOLOCAL) -np $(NP) $(MCH) $(PROG) -case $(CASE).depl

#w Makes local link for running with debugger
link:
	$(MAKE) -C $(PETSCFEM_DIR) $(LOC_PROG)
	ln -f $(PROG) $(LOC_PROG)
