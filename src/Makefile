#__INSERT_LICENSE__
# $Id: Makefile,v 1.40.2.7 2002/01/09 16:31:07 mstorti Exp $ 

.PHONY: all run lclean save compile

PETSCFEM_DIR = ..

# Add new files here
LIBOBJS_ = gpdata texthash readmesh elemset elemlist \
          fastmat utils getprop \
	  dofmap fem idmap arglist version util2 \
	  fstack tempfun fastmat2 fm2eperl \
	  measperf mainutl sttfilter state secant getarray \
	distmat sparse randomg spsolve pfmatsp texthf timestat \
	metisprt debug graph iisdcr iisdgraph iisdmat petscmat

# C files
CLIBOBJS_ =  getarrgr.tab

# Fortran files
FLIBOBJS_ =  spline seval

spline.o: spline.f
ssval.o: seval.f

LIBOBJS := $(LIBOBJS_:=.o) $(FLIBOBJS_:=.o) $(CLIBOBJS_:=.o) 

.INTERMEDIATE: $(LIBOBJS) 

SRCS := $(LIBOBJS_:=.cpp) $(CLIBOBJS_:=.c) 

TAGFLAGS = 

include $(PETSCFEM_DIR)/Makefile.base

all: libpetscfem.a

compile:  libpetscfem.a

libpetscfem.a: libpetscfem.a($(LIBOBJS))
	ranlib $@

libpetscfem.a(%.o): %.o
	ar ru $@ $<

#libpetscfem($(LIBOBJS)): $(LIBOBJS)
#libpetscfem.a: $(LIBOBJS)
#	ar -ru $@ $(LIBOBJS)
#	ranlib $@

#spline.o: spline.f
#	g77 -c $(FFLAGS) -O0 $<

.SECONDARY: fm2eperl.cpp

fm2eperl.o: fm2eperl.cpp
fm2eperl.cpp: fmat2ep.cpp prep.pl
	if [ -e $@ ] ; then chmod +w $@ ; rm $@ ; fi
	eperl -B '//<' -E '//>' $< > $@ ; chmod -w $@

FSM = matFSM
$(FSM).cpp $(FSM).h: sparse.sm
	$(SMC) < $<
	mv -f $(FSM).cc $(FSM).cpp

pfmat.sm: pfmat.sm.in
	if [ -f $@ ] ; then chmod +w $@; fi ;	\
	eperl -P $< > $@ ;			\
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
measperf.o: measperf.cpp
	-${CC} -c ${COPTFLAGS} ${CFLAGS} ${CCPPFLAGS} -O1 $<

local_clean:: local_dist_clean
	-rm -f version.cpp libpetscfem.a # &> /dev/null

readlist.h: readlist.eperl
	if [ -e readlist.h ] ; then chmod +w readlist.h ; fi
	eperl readlist.eperl > readlist.h
	chmod -w readlist.h

version.cpp: version.cppin ../VERSION
	rm -f version.cpp
	eperl version.cppin >version.cpp
	-chmod -w version.cpp

.SECONDARY: tryme.o 
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
	$(GNU_C_COMPILER) -c $(CPPFLAGS) $(CFLAGS) $< -o $@

odocsrc: $(DOC_DIR)/odociisd.tex $(DOC_DIR)/odocelems.tex $(DOC_DIR)/odocrmsh.tex

$(DOC_DIR)/odociisd.tex: lusubd.cpp
	$(run_odoc)

$(DOC_DIR)/odocelems.tex: elemset.cpp 
	$(run_odoc)

$(DOC_DIR)/odocrmsh.tex: readmesh.cpp
	$(run_odoc)
