#__INSERT_LICENSE__
 # $Id: Makefile,v 1.22 2001/05/23 20:23:42 mstorti Exp $ 

.PHONY: all run lclean save compile

PROG = advdif
PETSCFEM_DIR   = ../..
MYOBJS = advdife.o advdifeo.o advdif.o ffadvfm2_old.o  ffburg.o ffswfm2t.o \
		wallswt.o bccadvfm2.o bccadvfm2_old.o  ffadvfm2.o ffburg_old.o \
		enthalpy.o genload.o 

TAGFLAGS = $(TAGINCL)="../../src/TAGS"

SRCS := $(MYOBJS:.o=.cpp)

include $(PETSCFEM_DIR)/Makefile.base

compile: $(PROG).bin

ffeuler.o: ffeuler.cpp
	-${CC} -c ${COPTFLAGS} ${CFLAGS} ${CCPPFLAGS} -O0 $<


# documentation of options with the script `odoc.pl'
odocadvdif: ../../doc/odocadvdif.tex  ../../doc/odocadvdife.tex \
	../../doc/odocadvfm2.tex ../../doc/odocburgers.tex 

dd := ../../doc

$(dd)/odocadvdif.tex: advdif.cpp $(dd)/advdif.tex
	$(run_odoc)

../../doc/odocadvdife.tex: advdife.cpp  $(dd)/advdif.tex
	$(run_odoc)

../../doc/odocadvfm2.tex: ffadvfm2.cpp  $(dd)/advdif.tex $(dd)/advdifop.tex
	$(run_odoc)

../../doc/odocburgers.tex: ffburg.cpp  $(dd)/advdif.tex
	$(run_odoc)

../../doc/odocswfm2t.tex: ffswfm2t.cpp  $(dd)/advdif.tex
	$(run_odoc)

local_dist_clean::
	-rm -f *.o
