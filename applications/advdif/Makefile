#__INSERT_LICENSE__
 # $Id: Makefile,v 1.40.2.1 2003/10/16 19:07:15 mstorti Exp $ 

.PHONY: all run lclean save compile

PROG = advdif
PETSCFEM_MODULE = Advective-Diffusive
PETSCFEM_DIR   = ../..
MYOBJS = 

LIBOBJS = advdife.o advdifeo.o advdif.o ffadvfm2_old.o  ffburg.o ffswfm2t.o	\
		wallswt.o bccadvfm2.o bccadvfm2_old.o  ffadvfm2.o ffburg_old.o	\
		enthalpy.o genload.o linhff.o diff.o aquifer.o stream.o		\
		blessf.o bubbly.o advec.o compdate.o gasflow.o advdfgth.o	\
		advdif_bubbly.o advhookf.o smoke.o ctff.o streamsw1d.o id.o nonlres.o

.INTERMEDIATE: $(LIBOBJS) 

TAGFLAGS = $(TAGINCL)="../../src/TAGS"

SRCS := $(MYOBJS:.o=.cpp) $(LIBOBJS:.o=.cpp)

include $(PETSCFEM_DIR)/Makefile.base

# ffeuler.o: ffeuler.cpp
# 	-${CC} -c ${COPTFLAGS} ${CFLAGS} ${CCPPFLAGS} -O0 $<

advdif_bubbly.o : CCPPFLAGS += -O0

dd := ../../doc

# documentation of options with the script `odoc.pl'
odocadvdif: $(dd)/odocadvdif.tex $(dd)/odocadvdif2.tex $(dd)/odocadvdife.tex \
	$(dd)/odocadvfm2.tex $(dd)/odocburgers.tex \
	$(dd)/odocswfm2t.tex $(dd)/odocgenl.tex 

$(dd)/odocadvdif.tex: advdif.cpp $(dd)/advdif.tex
	$(run_odoc)

$(dd)/odocadvdif2.tex: streamsw1d.cpp $(dd)/advdif2.tex
	$(run_odoc)

$(dd)/odocadvdife.tex: advdife.cpp  $(dd)/advdif.tex
	$(run_odoc)

$(dd)/odocadvfm2.tex: ffadvfm2.cpp  $(dd)/advdif.tex $(dd)/advdifop.tex
	$(run_odoc)

$(dd)/odocburgers.tex: ffburg.cpp  $(dd)/advdif.tex
	$(run_odoc)

$(dd)/odocswfm2t.tex: ffswfm2t.cpp  $(dd)/advdif.tex
	$(run_odoc)

$(dd)/odocgenl.tex: genload.cpp linhff.cpp
	-chmod +w $@
	$(ODOC) -o $@ $^
	-chmod -w $@

$(dd)/odocdiff.tex: diff.cpp
	$(run_odoc)

$(dd)/odocstream.tex: stream.h stream.cpp
	-chmod +w $@
	$(ODOC) -o $@ $^
	-chmod -w $@

local_dist_clean::
	-rm -f *.o
