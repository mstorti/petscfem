# $Id: Makefile,v 1.9 2001/01/22 12:58:24 mstorti Exp $ 

.PHONY: all run lclean save compile

PETSCFEM_DIR = ..

# Add new files here
LIBOBJS_ = gpdata texthash readmesh elemset elemlist \
          fastmat utils getprop \
	  dofmap fem idmap arglist version util2 \
	  fstack tempfun fastmat2 fm2eperl \
	  measperf mainutl sttfilter state

# Fortran files
FLIBOBJS_ =  spline seval
spline.o: spline.f
ssval.o: seval.f

LIBOBJS := $(LIBOBJS_:=.o) $(FLIBOBJS_:=.o)

.INTERMEDIATE: $(LIBOBJS)

SRCS := $(LIBOBJS_:=.cpp)

TAGFLAGS = 

include $(PETSCFEM_DIR)/Makefile.base

all: libpetscfem.a

compile:  libpetscfem.a

libpetscfem.a: libpetscfem.a($(LIBOBJS))
	ranlib $@

libpetscfem.a(%.o): %.o
	ar -ru $@ $<

#libpetscfem($(LIBOBJS)): $(LIBOBJS)
#libpetscfem.a: $(LIBOBJS)
#	ar -ru $@ $(LIBOBJS)
#	ranlib $@

#spline.o: spline.f
#	g77 -c $(FFLAGS) -O0 $<

.SECONDARY: fm2eperl.cpp fm2eperl.o

fm2eperl.o: fm2eperl.cpp
fm2eperl.cpp: fmat2ep.cpp prep.pl
	if [ -e $@ ] ; then chmod +w $@ ; rm $@ ; fi
	eperl -B '//<' -E '//>' $< > $@ ; chmod -w $@

# This doesn't compile with -O2
measperf.o: CCPPFLAGS += -O1
#measperf.o: measperf.cpp
#	-${CC} -c ${COPTFLAGS} ${CFLAGS} ${CCPPFLAGS} -O1 $<

local_clean:: local_dist_clean
	-rm -f version.cpp libpetscfem.a # &> /dev/null

tryme: tryme.o seval.o spline.o
	g++ -o tryme tryme.o seval.o spline.o

readlist.h: readlist.eperl
	if [ -e readlist.h ] ; then chmod +w readlist.h ; fi
	eperl readlist.eperl > readlist.h
	chmod -w readlist.h

version.cpp: version.cppin ../VERSION
	rm -f version.cpp
	eperl version.cppin >version.cpp
	-chmod -w version.cpp
