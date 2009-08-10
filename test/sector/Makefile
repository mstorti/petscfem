#__INSERT_LICENSE__

.PHONY: constraint_bug chunk_size_process sector iprops triangle lap_per

default: sector

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include ../Makefile.base.test
SRCS = 

tests: sector chunk_size_process triangle lap_per

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
#w Test on an angular sector, uses constraints. 
sector: sector.depl $(LAPLACE)
	-$(MPIEXEC) -n 1 $(LAPLACE) \
		-case sector.depl $(call in_backg,sector)
	-mv sector.sal sector0.sal

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
NP := 2
MCH := -machinefile machi.dat 
#w Test chunk_size with the new chunk processing
chunk_size_process: $(LAPLACE)
	echo `uname -n` >machi.dat
	echo `uname -n` >>machi.dat
#
# 	Try with different chunk sizes and number of processors
#
	-chmod +w sector.depl 
	$(EPERL) -P -d chunk_size=1 sector.epl >sector.depl
	-$(MPIEXEC) $(MCH) -n $(NP) $(LAPLACE) -case sector.depl \
		$(call in_backg,$@-1)
	mv sector.sal sector_chunk1.sal
#
	-chmod +w sector.depl 
	$(EPERL) -P -d chunk_size=10 sector.epl >sector.depl
	-$(MPIEXEC) $(MCH) -n $(NP) $(LAPLACE) -case sector.depl \
		$(call in_backg,$@-2)
	mv sector.sal sector_chunk10.sal
#
	-chmod +w sector.depl 
	$(EPERL) -P -d chunk_size=100 sector.epl >sector.depl
	-$(MPIEXEC) $(MCH) -n $(NP) $(LAPLACE) -case sector.depl \
		$(call in_backg,$@-3)
	mv sector.sal sector_chunk100.sal
	-rm machi.dat

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
.PHONY: triangle triangle_npg3 triangle_npg4 triangle_npg7

#w Tests laplace with triangles
triangle: triangle_npg3 triangle_npg4 triangle_npg7

triangle_npg3:
	-chmod +w sector.depl 
	$(EPERL) -P -d triangle=1 -d npg=3 sector.epl >sector.depl
	-$(MPIEXEC) -n 1 $(LAPLACE) -case sector.depl \
		$(call in_backg,$@)
	mv sector.sal sector_triangle_npg3.sal

triangle_npg4:
	-chmod +w sector.depl 
	$(EPERL) -P -d triangle=1 -d npg=4 sector.epl >sector.depl
	-$(MPIEXEC) -n 1 $(LAPLACE) -case sector.depl  \
		$(call in_backg,$@)
	mv sector.sal sector_triangle_npg4.sal

triangle_npg7:
	-chmod +w sector.depl 
	$(EPERL) -P -d triangle=1 -d npg=7 sector.epl >sector.depl
	-$(MPIEXEC) -n 1 $(LAPLACE) -case sector.depl $(call in_backg,$@)
	mv sector.sal sector_triangle_npg7.sal

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
# checks integer properties
iprops:
	chmod +w sector.depl
	$(EPERL) -P -d chunk_size=10 -d iprops=1 sector.epl >sector.depl
	chmod -w sector.depl

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
#w Lapace on a square with periodic b.c.'s
lap_per: lap_per.dat laplace
	-$(MPIEXEC) -n 1 $(LAPLACE) -case lap_per.dat \
		$(call in_backg,$@)

LOC_PROG := laplace
PROG := $(LAPLACE)
CASE_NAME := sector

#
