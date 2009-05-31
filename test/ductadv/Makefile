# $Id cases-1.0.31-132-ga72f976 Thu Dec 6 16:33:22 2007 -0300$ 

.PHONY: all run  compile lclean laplace adv nso laplace_tests force \
        ctest gprof srfgath

default: tests

tests: checktmean

# BOPT := O_c++
PETSCFEM_DIR := ../../
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

LOC_PROG := advdif
CASE := ductadv
PROG := $(ADVDIF)

CASE := $(CASE)

################################################################
## Default
NP := 1
NL := 
MCH := 

checktmean: 
	mkstepsdir.pl
	make-depl $(CASE).epl $(CASE).depl -d mkmesh=1
	cp -f *.pl *.depl *.epl STEPS
	${MPI_HOME}/bin/mpiexec -a $(CASE) \
		$(MCH) -np $(NP) $(ADVDIF) \
		-case $(CASE).depl $(call redir,$@.output.tmp)
	$(MAKE) zip
	$(OCTAVE) checktmean.m > verif.tmp

prog:
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
	find STEPS* -name '$(CASE).state*.tmp' |	\
		xargs -r gzip -f 

#w Make link for debugging
link:
	$(MAKE) BOPT=g_c++ -C $(PETSCFEM_DIR) advdif
	ln -sf $(ADVDIF_G) .

kill:
	mpdkilljob -a $(CASE)
