#__INSERT_LICENSE__
# $Id merge-with-petsc-233-53-g55c6e64 Fri Oct 26 12:06:39 2007 -0300$ 

.PHONY: all run lclean save compile default

default: compile

PETSCFEM_DIR = ..

# Add new files here
LIBOBJS_ = \
amplidl arglist autostr blessf cloud cloud2 debug distmat dlhook	\
dofmap dofmap2 dofmap3 dvecpar dvector dxelmst dxhook elemlist		\
elemset elmsupl elmsupl2 fastmat fastmat2 fastmat2b fastmat2c		\
fastmat2d fem fm2ctx2 fm2eperl fm2eperl2 fm2temp fstack gatherer	\
generror getarray getprop gpdata graph graphdv hook idmap iisdcr	\
iisdgraph iisdmat iisdmat2 interpola lagmul linkgraph mainutl		\
matchg measperf metisprt penalize petscmat pfmat pfobject propef	\
randomg readmesh rhhook secant shllhook sparse spdirect spetscmat	\
spsolve srfgath state sttfilter surf2vol syncbuff tempfun		\
tempfun2 texthash texthf timestat util2 util3 utils fm2prod2

# C files
CLIBOBJS_ =  getarrgr.tab mkl_gf

# Fortran files
FLIBOBJS_ =  spline seval dgefa

GFORTRAN = gfortran -fPIC
spline.o: spline.f
	${GFORTRAN} -c spline.f
seval.o: seval.f
	${GFORTRAN} -c seval.f
dgefa.o: dgefa.f
	${GFORTRAN} -c dgefa.f

#spline.o: spline.f
#seval.o: seval.f
#dgefa.o: dgefa.f

LIBOBJS := $(LIBOBJS_:=$(osfx).o) $(FLIBOBJS_:=$(osfx).o) \
	$(CLIBOBJS_:=$(osfx).o) 

.INTERMEDIATE: $(LIBOBJS) 

SRCS := $(LIBOBJS_:=.cpp) $(CLIBOBJS_:=.c) 

TAGFLAGS = 

include $(PETSCFEM_DIR)/Makefile.base

LIBPETSCFEM := libpetscfem$(osfx).a

all: compile

compile:  $(LIBPETSCFEM)

$(LIBPETSCFEM): $(LIBPETSCFEM)($(LIBOBJS))
	ranlib $@

$(LIBPETSCFEM)(%.o): %.o
	ar ru $@ $<
	-rm -f $<

.SECONDARY: fm2eperl.cpp

fm2eperl$(osfx).o: fm2eperl.cpp
fm2eperl.cpp: fmat2ep.cpp prep.pl
	if [ -e $@ ] ; then chmod +w $@ ; rm $@ ; fi
	$(EPERL) -B '//<' -E '//>' $< > $@ ; chmod -w $@

fm2eperl2$(osfx).o: fm2eperl2.cpp
fm2eperl2.cpp: fmat2ep2.cpp prep.pl fmat2.pl
	if [ -e $@ ] ; then chmod +w $@ ; rm $@ ; fi
	$(EPERL) -B '//<' -E '//>' $< > $@ ; chmod -w $@

fem2eperl3$(osfx).o: fm2eperl3.cpp
fm2eperl3.cpp: fmat2ep3.cpp prep.pl fmat2.pl
	if [ -e $@ ] ; then chmod +w $@ ; rm $@ ; fi
	$(EPERL) -B '//<' -E '//>' $< > $@ ; chmod -w $@

## I disable this (enable if neccessary)
## generation of FSM files is left to `developper' mode
## ifeq ($(DEVELOPER_MODE),yes)
ifeq (1,0)
FSM := matFSM
$(FSM).cpp $(FSM).h: sparse.sm
	$(SMC) < $<
	mv -f $(FSM).cc $(FSM).cpp

pfmat.sm: pfmat.sm.in
	if [ -f $@ ] ; then chmod +w $@; fi ;	\
	$(EPERL) -P $< > $@ ;			\
	chmod -w $@

pfmatFSM.cpp  pfmatFSM.h: pfmat.sm
%FSM.cpp %FSM.h: %.sm
	if [ -f $*FSM.cpp ] ; then chmod +w $*FSM.cpp ; fi 
	if [ -f $*FSM.h ] ; then chmod +w $*FSM.h ; fi 
	$(SMC) < $<
	-rm -f $*FSM.cpp
	insdeb.pl $*FSM.cc $*FSM.cpp > /dev/null
	if [ ! -f $*FSM.cpp ] ; then cp -f $*FSM.cc $*FSM.cpp ; fi
	-rm -f $*FSM.cc
	chmod -w $*FSM.cpp $*FSM.h

readlist.h: readlist.eperl
	if [ -e readlist.h ] ; then chmod +w readlist.h ; fi
	$(EPERL) readlist.eperl > readlist.h
	chmod -w readlist.h

mprod.h mproddef.h: mprod.pl mprod1.h mprodtmpl.h
	./mprod.pl

endif

# This doesn't compile with -O2
#measperf.o: CCPPFLAGS += -O1 // this doesn't work with the old 'make'
# measperf.o: measperf.cpp
#	-${CC} -c ${COPTFLAGS} ${CFLAGS} ${CCPPFLAGS} -O1 $<

local_clean:: local_dist_clean
	-rm -f version.cpp $(LIBPETSCFEM) # &> /dev/null
pfmat.o dofmap2.o: CCPPFLAGS += -O0

version.cpp: version.cppin ../VERSION
	rm -f version.cpp
	$(EPERL) version.cppin >version.cpp
	-chmod -w version.cpp

version.o: CPPFLAGS += -DPETSC_BOPT_FLAG="\"$(BOPT)\""

#.SECONDARY: tryme$(osfx).o 
#tryme: CPPFLAGS = 
#tryme: secant.o tryme.o
#	g++ -o tryme tryme.o secant.o 
#secant.o: secant.cpp secant.h
#tryme.o: tryme.cpp secant.h

#secant.o: secant.cpp
#	g++ -fPIC -c $<

#tryme.o: tryme.cpp
#	g++ -c $<

## `getarrgr.y' is rarely modified, it's pretty frozen.
## In case it is modified we have to uncomment the 
## following line. 
## ifeq ($(DEVELOPER_MODE),yes)
ifeq (1,0)
local_clean::
	-rm -rf tryme *.tab.*

local_depend: getarrgr.tab.c getarrgr.tab.h

#.INTERMEDIATE: *.tab.*
%.tab.c %.tab.h: %.y
	$(BISON) $<

getarrgr.tab.o: getarrgr.tab.c
endif

odocsrc: $(DOC_DIR)/odociisd.tex $(DOC_DIR)/odociisdm.tex \
		$(DOC_DIR)/odocelems.tex $(DOC_DIR)/odocrmsh.tex \
		$(DOC_DIR)/odocdxhook.tex

$(DOC_DIR)/odociisd.tex: iisdcr.cpp
	$(run_odoc)

$(DOC_DIR)/odociisdm.tex: iisdmat.cpp
	$(run_odoc)

$(DOC_DIR)/odocelems.tex: elemset.cpp 
	$(run_odoc)

$(DOC_DIR)/odocrmsh.tex: readmesh.cpp
	$(run_odoc)

$(DOC_DIR)/odocdxhook.tex: dxhook.cpp
	$(run_odoc)

PFVH := $(PETSCFEM_DIR)/tools/pfversion.h.undef
pfversion.h: force
	if ! docdate.pl -l C -o pfversion.h ; then cp $(PFVH) pfversion.h ; fi
