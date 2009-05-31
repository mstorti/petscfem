# $Id cases-1.0.31-132-ga72f976 Thu Dec 6 16:33:22 2007 -0300$ 

.PHONY: all run  compile lclean laplace adv nso laplace_tests force \
        ctest gprof srfgath

default: runadv

BOPT := O_c++
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

LOC_PROG := ns
CASE := ductadv
PROG := $(NS)

CASE := $(CASE)

################################################################
## Default
NP := 1
NL := 
MCH := 

run: 
	mkstepsdir.pl
	make-depl $(CASE).epl $(CASE).depl -d mkmesh=1
	cp -f *.pl *.depl *.epl STEPS
ifeq (1,1)
	${MPI_HOME}/bin/mpiexec -a $(CASE) \
		$(MCH) -np $(NP) $(ADVDIF) -case $(CASE).depl
	$(MAKE) zip
endif

cylinder.nod-vel.tmp: cylinder.state.tmp

cpydata:

prog:
#	$(MAKE) -C $(PETSCFEM_DIR) ns
	$(MAKE) -C $(PETSCFEM_DIR) advdif

local_clean::
	-chmod +w *.depl
	-rm *.depl *.sal tmp_file*.tmp &> /dev/null
	shopt -s nullglob ; rrm -rf STEPS*

shell_init:

shell_pre:

shell_post: 
	shopt -s nullglob ;				\
	find STEPS* -name 'cylinder*.state*.tmp' |	\
		xargs -r gzip -f 

shell_close:

zip:
	shopt -s nullglob ;				\
	find STEPS* -name 'cylinder*.state*.tmp' |	\
		xargs gzip -f 

#w Make link for debugging
link:
	$(MAKE) BOPT=g_c++ -C $(PETSCFEM_DIR) advdif
	ln -sf $(ADVDIF_G) .

kill:
	mpdkilljob -a $(CASE)

DIRE := /u/mstorti/PETSC/petscfem-cases/animation/particle-hips/cylinder-rtd
stcpy:
	scp coyotep:$(DIRE)/STEPS/cylinder-adv.state-200.tmp.gz cylinder-adv.state.tmp.gz 
	gunzip -f cylinder-adv.state.tmp.gz 
