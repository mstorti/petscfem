# $Id: Makefile,v 1.1.1.1 2000/12/28 12:54:43 mstorti Exp $ 

.PHONY: constraint_bug chunk_size_process sector iprops

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include ../Makefile.base.test
SRCS = 

# with this it doesn't work
ADVO    = $(PETSCFEM_DIR)/../petscfem-1.48g/applications/advective/adv.bin 

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
sector: sector.depl $(LAPLACE)
	-$(MPIRUN) -np 1 $(LAPLACE) -case sector.depl  >/dev/null
	-mv sector.sal sector0.sal

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
# Test chunk_size with the new chunk processing
NP := 2
MCH := -machinefile machi.dat 
chunk_size_process: $(LAPLACE)
	echo `uname -n` >machi.dat
	echo `uname -n` >>machi.dat
#
# 	Try with different chunk sizes and number of processors
#
	-chmod +w sector.depl 
	eperl -P -d chunk_size=1 sector.epl >sector.depl
	-$(MPIRUN) $(MCH) -np $(NP) $(LAPLACE) -case sector.depl >/dev/null
	mv sector.sal sector_chunk1.sal
#
	-chmod +w sector.depl 
	eperl -P -d chunk_size=10 sector.epl >sector.depl
	-$(MPIRUN) $(MCH) -np $(NP) $(LAPLACE) -case sector.depl >/dev/null
	mv sector.sal sector_chunk10.sal
#
	-chmod +w sector.depl 
	eperl -P -d chunk_size=100 sector.epl >sector.depl
	-$(MPIRUN) $(MCH) -np $(NP) $(LAPLACE) -case sector.depl >/dev/null
	mv sector.sal sector_chunk100.sal
	-rm machi.dat

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
.PHONY: triangle triangle_npg3 triangle_npg4 triangle_npg7

triangle: triangle_npg3 triangle_npg4 triangle_npg7

triangle_npg3:
	-chmod +w sector.depl 
	eperl -P -d triangle=1 -d npg=3 sector.epl >sector.depl
	-$(MPIRUN) $(MCH) -np 1 $(LAPLACE) -case sector.depl  >/dev/null
	mv sector.sal sector_triangle_npg3.sal

triangle_npg4:
	-chmod +w sector.depl 
	eperl -P -d triangle=1 -d npg=4 sector.epl >sector.depl
	-$(MPIRUN) $(MCH) -np 1 $(LAPLACE) -case sector.depl  >/dev/null
	mv sector.sal sector_triangle_npg4.sal

triangle_npg7:
	-chmod +w sector.depl 
	eperl -P -d triangle=1 -d npg=7 sector.epl >sector.depl
	-$(MPIRUN) $(MCH) -np 1 $(LAPLACE) -case sector.depl  >/dev/null
	mv sector.sal sector_triangle_npg7.sal

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
# checks integer properties
iprops:
	chmod +w sector.depl
	eperl -P -d chunk_size=10 -d iprops=1 sector.epl >sector.depl
	chmod -w sector.depl

include makefile.d
