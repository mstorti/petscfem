# $Id: Makefile,v 1.2 2001/01/02 22:28:18 mstorti Exp $ 

.PHONY: homof sqcav21 force
run: sqcav

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
NS = $(PETSCFEM_DIR)/applications/ns/ns.bin

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
TMP_FILES = homof.con.tmp  homof.ini.tmp  homof.nod.tmp  homof.peri.tmp 

homof: homof.depl $(TMP_FILES)
	-$(MPIRUN) -np 1 $(NS) -case homof.depl # >/dev/null

$(TMP_FILES): mkmesh.m
	octave -q < $<

.PHONY: lclean add_clean

local_clean::
	chmod +w *.depl
	-@rm output.m a.dat *.sal *.depl &>/dev/null

CASE = back_step

ns:
	make -C $(PETSCFEM_DIR) ns

sqcav21: sqcav21.depl ns
	-$(MPIRUN) -np 1 $(NS) -case sqcav21.depl  | tee salida.sal

MCH = -machinefile ./machi.dat
sqcav: sqcav.depl ns
	-procsel proctable
	-$(MPIRUN) -np 1 $(MCH) $(NS) -case sqcav.depl  | tee salida.sal

$(CASE).state: 
	-chmod +w $(CASE).depl
	eperl -P -d verify=0 $<   > $(CASE).depl
	-chmod -w $(CASE).depl
	-$(MPIRUN) -np 1 $(NS) -case $(CASE).depl  > salida.sal

system.dat: $(CASE).state
	-chmod +w $(CASE).depl
	eperl -P -d verify=1 $(CASE).epl   > $(CASE).depl
	-chmod -w $(CASE).depl
	-$(MPIRUN) -np 1 $(NS) -case $(CASE).depl  >salida.sal

.PHONY: dat
dat: mat.dat id.dat

mat.dat: salida.sal system.dat
	readdata.pl

