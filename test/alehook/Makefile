# $Id: Makefile,v 1.1 2005/01/07 01:44:20 mstorti Exp $ 

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test
ALEHOOK := "~/PETSC/ALEHOOK/alehook.efn";

LOC_PROG := ns
CASE_NAME := sqcav
PROG := $(NS)


###--------<*>--------<*>--------<*>--------<*>--------<*>--------<*>--------<*> 
NP = 1
MCH = -machinefile ./machi.dat

#w Runs program
run:
	-${MPIRUN} -np $(np) $(MCH) $(PROG) -case $(CASE).depl \
	 	      	   > $(output_file)

