# $Id: Makefile,v 1.2 2005/03/28 16:42:59 mstorti Exp $ 

default: condwall

# BOPT = g_c++
PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

LOC_PROG := ns
PROG := $(NS)

NP = 1
# MCH = -machinefile ./machi.dat
# NP := $(shell procsel $(PETSCFEM_PROCTABLE))

.PHONY: default condwall

condwall:
	$(MAKE) CASE=condwall run

#w Runs program
run:
	make-depl $(CASE).epl $(CASE).depl
	-${MPIRUN} $(MCH) -np $(NP) $(PROG) -case $(CASE).depl

###--------<*>--------<*>--------<*>--------<*>--------<*>--------<*>
#w Makes local link for running with debugger
link:
	$(MAKE) -C $(PETSCFEM_DIR) $(LOC_PROG)
	ln -f $(PROG) $(LOC_PROG)

