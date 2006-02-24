# $Id: Makefile,v 1.1.30.1 2006/02/24 12:26:57 rodrigop Exp $ 

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test
ALEHOOK := "~/PETSC/ALEHOOK/alehook.efn";

LOC_PROG := ns
CASE_NAME := sqcav
PROG := $(NS)


###--------<*>--------<*>--------<*>--------<*>--------<*>--------<*>--------<*> 
NP = 1
MCH =    

#w Runs program
run:
	-${MPIRUN} -np $(np) $(MCH) $(PROG) -case $(CASE).depl \
	 	      	   > $(output_file)

