# $Id: Makefile,v 1.1.1.1 2000/12/28 12:54:42 mstorti Exp $ 

.PHONY: all run lclean save compile

PETSCFEM_DIR = ..
include $(PETSCFEM_DIR)/Makefile.base

LIBOBJS = gpdata.o texthash.o readmesh.o elemset.o elemlist.o \
          fastmat.o utils.o getprop.o \
	  dofmap.o fem.o idmap.o arglist.o version.o util2.o \
	  fstack.o tempfun.o fastmat2.o fm2eperl.o spline.o seval.o \
	  measperf.o


SRCS := $(LIBOBJS:.o=.cpp)

TAGFLAGS = 

all: libpetscfem.a

compile:  libpetscfem.a

libpetscfem.a: $(LIBOBJS)
	ar -ru $@ $(LIBOBJS)
	ranlib $@

#spline.o: spline.f
#	g77 -c $(FFLAGS) -O0 $<

.SECONDARY: fm2eperl.cpp fm2eperl.o

fm2eperl.o: fm2eperl.cpp
fm2eperl.cpp: fmat2ep.cpp prep.pl
	if [ -e $@ ] ; then chmod +w $@ ; rm $@ ; fi
	eperl -B '//<' -E '//>' $< > $@ ; chmod -w $@

# This doesn't compile with -O2
measperf.o: measperf.cpp
	-${CC} -c ${COPTFLAGS} ${CFLAGS} ${CCPPFLAGS} -O1 $<


lclean: general_clean
	-rm *.o libpetscfem.a &> /dev/null
#
# Esto deberia ir en realclean o algo por el estilo
#	if [ -e fm2eperl.cpp ] ; then chmod +w fm2eperl.cpp ; \
#		rm fm2eperl.cpp ; fi

tag:
	@cd $(HOME)/PETSC/petscfem ; cvs history -T ; \
	echo -n "enter tag: " ; read TAG ; \
	echo -n "enter comment: " ; read COMMENT ; \
	echo "setting tag $$TAG ..." ; \
	cvs rtag $$TAG -m "$$COMMENT" petscfem ; \
	echo " ...  done."

tryme: tryme.o seval.o spline.o
	g++ -o tryme tryme.o seval.o spline.o

readlist.h: readlist.eperl
	if [ -e readlist.h ] ; then chmod +w readlist.h ; fi
	eperl readlist.eperl > readlist.h
	chmod -w readlist.h

version.cpp: version.cppin $(TARFILE).tara,v
	fversion -il ../save.log $(TARFILE).tara,v \
			version.cppin version.cpp

# File containing dependencies computed by makedepend
include makefile.d
