# $Id: Makefile,v 1.4 2005/03/29 04:01:36 mstorti Exp $ 

default: condwall

# BOPT = g_c++
PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

LOC_PROG := ns
PROG := $(NS)
CASE := condwall

NP = 1
# MCH = -machinefile ./machi.dat
# NP := $(shell procsel $(PETSCFEM_PROCTABLE))

.PHONY: default condwall

condwall:
	$(MAKE) run

###--------<*>--------<*>--------<*>--------<*>--------<*>--------<*>
#w Makes local link for running with debugger
movwall.efn: movwall.cpp $(LIBPETSCFEM)

#w Runs program
run: movwall.efn
	make-depl $(CASE).epl $(CASE).depl
	-${MPIRUN} $(MCH) -np $(NP) $(PROG) -case $(CASE).depl

###--------<*>--------<*>--------<*>--------<*>--------<*>--------<*>
#w Runs PETSc-FEM instance for DX
dx:
	make-depl $(CASE).epl $(CASE)-dx.depl -d dx=1
	-${MPIRUN} -np 1 $(PROG) -case $(CASE)-dx.depl

zip:
	for f in STEPS/$(CASE).state_*.tmp ; \
	do echo gzipping $$f ; gzip -f $$f ; done

dx_make_command:
	octave -qH fixstate.m

###--------<*>--------<*>--------<*>--------<*>--------<*>--------<*>
#w Makes local link for running with debugger
link:
	$(MAKE) -C $(PETSCFEM_DIR) $(LOC_PROG)
	ln -f $(PROG) $(LOC_PROG)

