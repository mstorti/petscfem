#__INSERT_LICENSE__
 # $Id: Makefile,v 1.59 2006/12/24 03:10:22 mstorti Exp $ 

.PHONY: all run lclean save compile

PROG = advdif
PETSCFEM_MODULE = Advective-Diffusive
PETSCFEM_DIR   = ../..

PROGOBJ = main.o

LIBOBJS = advdif.o advdif_bubbly.o advdif_fsi.o compdate.o		\
	advdifglb.o advdife.o advdifeo.o  ffadvfm2_old.o ffburg.o	\
	ffswfm2t.o wallswt.o ffadvfm2.o	\
	ffburg_old.o enthalpy.o genload.o linhff.o diff.o aquifer.o	\
	stream.o blessf.o bubbly_full.o bubbly_asm.o advec.o		\
	gasflow.o advdfgth.o  advhookf.o smoke.o ctff.o			\
	streamsw1d.o id.o gaschem.o nonlres.o fm2funm.o			\
	bubbly_asm_mix.o lagmul.o advabso.o gasflwgth.o svenant1d.o	\
	advpenal.o absolay.o advdife-gcl.o advdife-ale.o streamsw2dt.o  \
	streamsw2d.o strat_sw2d.o advdife_dual_time.o gasflow_preco.o	\
	streamsw2d.o strat_sw2d.o advabsow.o advdif_dual_time.o         \
	bccadvfm2.o bccadvfm2_old.o volabso.o advdife-bdf.o advdif-bdf.o

MYOBJS = 

.INTERMEDIATE: $(LIBOBJS) 

TAGFLAGS = $(TAGINCL)="../../src/TAGS"

SRCS := $(MYOBJS:.o=.cpp) $(LIBOBJS:.o=.cpp)

include $(PETSCFEM_DIR)/Makefile.base

# ffeuler.o: ffeuler.cpp
# 	-${CC} -c ${COPTFLAGS} ${CFLAGS} ${CCPPFLAGS} -O0 $<

advdif_bubbly.o : CCPPFLAGS += -O0

dd := ../../doc

# documentation of options with the script `odoc.pl'
odocadvdif: $(dd)/odocadvdif.tex		\
	$(dd)/odocadvdif2.tex			\
	$(dd)/odocadvdife.tex			\
	$(dd)/odocadvfm2.tex			\
	$(dd)/odocburgers.tex			\
	$(dd)/odocswfm2t.tex			\
	$(dd)/odocgenl.tex			\
	$(dd)/odocstream.tex			\
	$(dd)/odocabso.tex  			\
	$(dd)/odocabsow.tex  			\
	$(dd)/odocdiff.tex  

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

$(dd)/odocabso.tex: advabso.cpp
	$(run_odoc)

$(dd)/odocabsow.tex: advabsow.cpp
	$(run_odoc)

local_dist_clean::
	-rm -f *.o
