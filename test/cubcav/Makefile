# $Id: Makefile,v 1.106 2004/11/09 19:41:51 mstorti Exp $ 

.PHONY: all run  compile lclean laplace adv nso laplace_tests force \
        ctest gprof srfgath

default: test_perfo

#default case
# case := laplace
# case := iisd_sbp1_uj
export case

tests: block_uploading

BOPT := O_c++
PETSCFEM_DIR = ../..
# PETSCFEM_DIR = $(HOME)/PETSC/petscfem-stable
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

LOC_PROG := ns
#CASE_NAME := qharm
# CASE_NAME := cubcav
CASE_NAME := cubcav-fs
PROG := $(NS)

CASE := $(CASE_NAME)
################################################################
ifeq ($(HOSTNAME),node1.beowulf.gtm)
NP := $(shell procsel -n ./proctable.geronimo)
NL := -nolocal
MCH := -machinefile ./machi.dat
else
NP := $(shell procsel ./proctable.mono)
NL :=
MCH := -machinefile ./machi.dat
NP := 1
endif

run_subdir_test:							\
	if [ -e $(sub_dir) ] ; then rrm -rf $(sub_dir) ; fi ;		\
	mkdir $(sub_dir) ;						\
	echo "node1 1.3 server" > $(sub_dir)/proctable ;		\
	for n in $(node_list) ;						\
		do grep $$n ./proctable >> $(sub_dir)/proctable ;	\
	done ;								\
	cp $(CASE_NAME).epl $(sub_dir) ;				\
	cp Makefile.sub-dir  $(sub_dir)/Makefile ;			\
	echo $(MAKE) -C $(sub_dir) data_dir=.. mkmesh=0

#w Makes only the mesh
mesh:
	$(MAKE) -C $(PETSCFEM_DIR)/tools hexasplit.bin
	make-depl $(CASE_NAME).epl $(CASE_NAME).depl -d case=ns -p '>>'

#w Runs program, use `mkmesh=0' or 1 to control mesh creation
run: 
	$(MAKE) mesh
	-${MPIRUN} $(NL) -np $(NP) $(MCH) $(PROG) -case $(CASE).depl

local_sw::
	chmod 755 mkmpeg*.sh

local_clean::
	-chmod +w *.depl
	-rm *.depl *.sal tmp_file*.tmp &> /dev/null

nettest = $(PETSCFEM_DIR)/tools/nettest.bin

.PHONY: nettest
nettest:
	NP=`procsel -n ./proctable.geronimo` ;  \
	$(MPIRUN) -np $$NP -machinefile machi.dat $(nettest)

dxo: viewer.dx.tmp

viewer.dx.tmp: viewer.dx
	$(EPERL) -P viewer.dx > viewer.dx.tmp

dx:
	$(EPERL) -P -d dx=1 $(CASE_NAME).epl > $(CASE_NAME)_dx.depl
	${MPIRUN} -np 1 $(PROG) -case $(CASE_NAME)_dx.depl

#================================================================
#w 
dx_make_command:
# dx_make_command: # dx_make_command2

dx_make_command1:
	octave -qH comp_step.m

dx_make_command2:
#	echo $$dx_step >> comp_server.fifo 
	if [ -f $(CASE_NAME).tiff ] ; 			\
		then rm $(CASE_NAME).tiff ; fi ;	\
	sendstep.pl $$dx_step				\
		comp_server.fifo comp_server2.fifo

#================================================================
yuv_compress2:

yuv_compress:
	yuv=YUV/$(CASE_NAME).$$dx_step.yuv ;	\
	tiff=$(CASE_NAME).tiff ;		\
	if [ -f $$tiff ] ; then			\
		if [ -f $$yuv ] ;		\
			then rm -f $$yuv ;	\
		fi ;				\
		convert $$tiff $$yuv ;		\
		gzip -f $$yuv ;			\
	fi

#================================================================
k:
	-killall ns_O.bin job.pl

#================================================================
#w Basic test for block uploading. This somewhat trivial since
#w for NS the profile is full. 
block_uploading:
	$(MAKE) CASE=$(CASE_NAME) case=plain_bupl0 run \
		> $(CASE_NAME).out.plain_bupl0.tmp
	$(MAKE) CASE=$(CASE_NAME) case=plain_bupl1 run \
		> $(CASE_NAME).out.plain_bupl1.tmp
	$(OCTAVE) check_plain.m				\
		> $(CASE_NAME).plain-bupl.verif.tmp

#---<*>---//---<*>---//---<*>---//---<*>---//---<*>---// 
#w Checks the surf-gatherer class
srfgath:
	$(MAKE) case=srfgath mesh
	-${MPIRUN} $(NL) -np $(NP) $(MCH) $(PROG) -case $(CASE).depl

#================================================================
#w For the shell-hook
cubcav_init:
	mkstepsdir.pl

cubcav_pre:
cubcav_post: zipall

zipall:
	for f in `find STEPS -name '*.state_*.tmp'` ; \
		do echo "gzipping $$f" ; gzip -f $$f ; done

cubcav_close:
#w
fifo:
	fifo="comp_server.fifo" ;			\
	if [ -e $$fifo ] ; then rm $$fifo ; fi ;	\
	mkfifo $$fifo
	fifo="comp_server2.fifo" ;			\
	if [ -e $$fifo ] ; then rm $$fifo ; fi ;	\
	mkfifo $$fifo

clean_fifo:
	fifo="comp_server.fifo" ;			\
	if [ -e $$fifo ] ; then rm $$fifo ; fi ;	\
	fifo="comp_server2.fifo" ;			\
	if [ -e $$fifo ] ; then rm $$fifo ; fi ;	\

#w Launch 
step_server:
	octave -qH comp_step2.m

#w Make link for debugging
link:
	$(MAKE) -C ../../ $(LOC_PROG)
	ln -f $(PROG) $(LOC_PROG)

#w performance test:
test_perfo:
	$(MAKE) N=132 run

sync:
	chmod 755 *.pl *.sh
