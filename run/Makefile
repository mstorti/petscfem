# $Id: Makefile,v 1.2 2001/01/02 22:28:14 mstorti Exp $ 

.PHONY: all run  compile lclean ns laplace adv nso laplace_tests force \
        ctest gprof

all: run

gprof: petscfem.prof

PETSCFEM_DIR = ..
include $(PETSCFEM_DIR)/Makefile.base

#PETSC_OPS = -trdump -log_info -trmalloc 
#PETSC_OPS = -log_info 
#PETSC_OPS = 

LAPLACE = $(PETSCFEM_DIR)/applications/laplace/laplace.bin
NS    = $(PETSCFEM_DIR)/applications/ns/ns.bin
ADV    = $(PETSCFEM_DIR)/applications/advective/adv.bin

#shortcuts
ns: $(NS)
adv: $(ADV)
laplace: $(LAPLACE)

$(NS): force
	$(MAKE) -C $(PETSCFEM_DIR) ns

$(ADV): force
	$(MAKE) -C $(PETSCFEM_DIR) adv

$(LAPLACE): force
	$(MAKE) -C $(PETSCFEM_DIR) laplace

################################################################
NP = 1
PROG  = $(NS)
DATAFILE = oscplate2.dat
MCH = -machinefile ./machi.dat
run: $(DATAFILE) force
	procsel proctable
	$(MAKE) $(PROG)
	-${MPIRUN} -np $(NP) $(MCH) $(PROG) -case $(DATAFILE) #&> salida.sal
# 	if [ -e tmp_file*.tmp ] ; then rm tmp_file*.tmp ; fi
# 	if [ -e gmon.out ] ; then $(MAKE) gprof ; fi

petscfem.prof: gmon.out
	@if [ -e gmon.out ] ; then echo "running gprof..." ; \
		gprof $(PROG) gmon.out >$@ ; fi

ctest: $(DATAFILE) force
	$(MAKE) $(NS)
	-${MPIRUN} -np 1 $(NS) -case oscplate.depl

tests: 
	-rm *.sal
	@echo "Running tests..."
	-${MPIRUN} -np 1 $(LAPLACE) -case lap_per.dat >/dev/null
	-${MPIRUN} -np 1 $(LAPLACE) -case sector.dat >/dev/null
	-${MPIRUN} -np 1 $(NS) -case oscplate.depl >/dev/null
	@echo "Verifying results..."
	runtests.pl
	@echo "Done."

naca.tmp.ini:
	perl -e 'require "../test/eperlini.pl"; \
		makeini([1.,1.,0.,180.5],2429,"$@");'

local_clean::
	-rm *.depl *.sal tmp_file*.tmp &> /dev/null
