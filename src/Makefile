#__INSERT_LICENSE__
# $Id: Makefile,v 1.92 2004/01/26 20:22:34 mstorti Exp $ 

.PHONY: all run lclean save compile default

default: compile

PETSCFEM_DIR = ..

# Add new files here
LIBOBJS_ = amplidl arglist debug distmat dofmap dofmap2 elemlist elemset	\
fastmat fastmat2 fem fm2eperl fm2eperl2 fstack getarray getprop gpdata		\
graph graphdv idmap iisdcr iisdgraph iisdmat linkgraph mainutl			\
measperf metisprt petscmat pfmat propef randomg readmesh secant sparse		\
spdirect spsolve state sttfilter tempfun tempfun2 texthash texthf		\
timestat util2 util3 utils version gatherer hook dlhook dxhook dxelmst autostr	\
cloud cloud2 surf2vol pfobject fm2temp shllhook dvector interpola elmsupl	\
dofmap3 iisdmat2 generror syncbuff srfgath blessf

# C files
CLIBOBJS_ =  getarrgr.tab

# Fortran files
FLIBOBJS_ =  spline seval dgefa

spline.o: spline.f
ssval.o: seval.f

LIBOBJS := $(LIBOBJS_:=$(osfx).o) $(FLIBOBJS_:=$(osfx).o) $(CLIBOBJS_:=$(osfx).o) 

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

FSM = matFSM
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

# This doesn't compile with -O2
#measperf.o: CCPPFLAGS += -O1 // this doesn't work with the old 'make'
# measperf.o: measperf.cpp
#	-${CC} -c ${COPTFLAGS} ${CFLAGS} ${CCPPFLAGS} -O1 $<

local_clean:: local_dist_clean
	-rm -f version.cpp $(LIBPETSCFEM) # &> /dev/null
readlist.h: readlist.eperl
	if [ -e readlist.h ] ; then chmod +w readlist.h ; fi
	$(EPERL) readlist.eperl > readlist.h
	chmod -w readlist.h
pfmat.o dofmap2.o: CCPPFLAGS += -O0


version.cpp: version.cppin ../VERSION
	rm -f version.cpp
	$(EPERL) version.cppin >version.cpp
	-chmod -w version.cpp
version.o: CPPFLAGS += -DPETSC_BOPT_FLAG="\"$(BOPT)\""


.SECONDARY: tryme$(osfx).o 
tryme: CPPFLAGS = 
tryme: secant.o tryme.o
	g++ -o tryme tryme.o secant.o 
secant.o: secant.cpp secant.h
tryme.o: tryme.cpp secant.h

secant.o: secant.cpp
	g++ -c $<

tryme.o: tryme.cpp
	g++ -c $<

local_clean::
	-rm -rf tryme *.tab.*

local_depend: getarrgr.tab.c getarrgr.tab.h

#.INTERMEDIATE: *.tab.*
%.tab.c %.tab.h: %.y
	$(BISON) $<

getarrgr.tab.o: getarrgr.tab.c

%.o: %.c
	$(GNU_C_COMPILER) -o $@ -c $(CPPFLAGS) $(CFLAGS) $< -o $@

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

pp:
	echo $(LIBPETSCFEM)
