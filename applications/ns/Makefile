#__INSERT_LICENSE__
 # $Id merge-with-petsc-233-52-g3182b2b Thu Oct 25 22:14:27 2007 -0300$ 

.PHONY: all run lclean save compile

PROG = ns
PETSCFEM_DIR   = ../..
PETSCFEM_MODULE = Navier-Stokes

PROGOBJ = main.o

LIBOBJS = ns.o mmvmain.o mmv2main.o nsstruct.o ns_fsi.o compdate.o	\
	nsglb.o nsitetlesfm2.o bccnsfm2.o wall.o nsvol.o walldata.o	\
	nsilesther.o bccnstherfm2.o nsikeps.o nsikepsrot.o		\
	wallke.o nonlr.o lagmul.o nssup.o nssupg.o nssupr.o nsid.o	\
	nsres.o adaptor.o elast.o adappg.o qharm.o qharmm.o		\
	nsirot.o blessf.o fracstep.o fstepfm2.o embgath.o		\
	mmove.o fm2funm.o genload.o nsgath.o nshookf.o invcoupl.o	\
	nullvort.o nsilesasm_2.o bccnsasmfm2.o bccfstep.o nslagmul.o	\
	condwall.o condwallpen.o nspenal.o nsroblgml.o fstepfm2cw.o	\
	nsitetlesc.o condwall-new.o mmove2.o nsitetasm.o		\
	bccnsitetasm.o nsitetasm_avgvol.o ns_gasflow.o bccns_gasflow.o	\
	mmoveopt.o elast2.o  bccnsitetasm_avgvol.o elastld.o		\
	errestim.o mmoveopt3.o bubblyqint.o nsitetlesf.o nsitetlesfbf.o	\
	poiboltz.o mmoveopt2.o truss.o nsitetlesfd.o nodeload.o		\
	nsitetlesls.o electrophoresis.o  electrophoresisM.o flowrev.o	\
	renorm.o renorm2.o renorm3.o nsitlfctx2.o					\
	electrophoresisM2.o electrophoresisM.o electrophoresis_mov.o	\
        poisson.o charge_cons.o pot_grad.o

MYOBJS =


TAGFLAGS = $(TAGINCL)="../../src/TAGS"
#TAGFLAGS = --include="../../src/TAGS"

# for makedepend
SRCS := $(MYOBJS:.o=.cpp) $(LIBOBJS:.o=.cpp)

include $(PETSCFEM_DIR)/Makefile.base

# CPPFLAGS += -DROSI_COUPLING_MODULE

# documentation of options with the script `odoc.pl'
odocns: ../../doc/odocns.tex ../../doc/odocnse.tex \
	../../doc/odocnsb.tex  ../../doc/odocnsw.tex \
	../../doc/odocnsfs.tex ../../doc/odocflrv.tex 

../../doc/odocns.tex: ns.cpp 
	$(run_odoc)

../../doc/odocnse.tex: nsitetlesfm2.cpp 
	$(run_odoc)

../../doc/odocnsb.tex: bccnsfm2.cpp 
	$(run_odoc)

../../doc/odocnsw.tex: wall.cpp 
	$(run_odoc)

../../doc/odocnsfs.tex: nssupg.cpp ../../doc/nsdoc.tex
	$(run_odoc)

../../doc/odocflrv.tex: flowrev.cpp 
	$(run_odoc)

ns.o : CCPPFLAGS += -O0

#compdate.o: CPPFLAGS += "-DPETSCFEM_HOSTNAME=\"$(shell uname -a | sed s/\#//g)\""
noafs:
	@echo "making symbolic links to NO Algebraic FS formulation... "
	@ls -l ns.cpp nsitetlesfm2.cpp ; \
	echo -n "re-link to .noafs files? (y/n) > " ; \
	read answer ; \
	if [ -L ns.cpp -a -L nsitetlesfm2.cpp -a "$$answer" = "y" ] ; \
		then ln -sf ns.cpp.noafs ns.cpp ; \
		ln -sf nsitetlesfm2.cpp.noafs nsitetlesfm2.cpp ; \
	fi ; 
	-@rm ns.o nsitetlesfm2.o

afs:
	@echo "making symbolic links to Algebraic FS formulation... "
	@ls -l ns.cpp nsitetlesfm2.cpp ; \
	echo -n "re-link to .afs files? (y/n) > " ; \
	read answer ; \
	if [ -L ns.cpp -a -L nsitetlesfm2.cpp -a "$$answer" = "y" ] ; \
		then ln -sf ns.cpp.afs ns.cpp ; \
		ln -sf nsitetlesfm2.cpp.afs nsitetlesfm2.cpp ; \
	fi ; 
	-@rm ns.o nsitetlesfm2.o

#lclean: general_clean
#	cat /dev/null >makefile.d
#	-rm *.o $(OBJS) $(PROG).bin octave-core *.cppi outvector.sal makefile.d
#	$(MAKE) clean

# File containing dependencies computed by makedepend
#
