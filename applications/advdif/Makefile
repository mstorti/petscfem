 # $Id: Makefile,v 1.11 2001/02/06 17:14:31 mstorti Exp $ 

.PHONY: all run lclean save compile

PROG = advdif
PETSCFEM_DIR   = ../..
MYOBJS = advdife.o advdif.o ffadvfm2.o  bccadvfm2.o  ffburg.o ffswfm2t.o \
		wall.cpp
#		absorb.o \
# 		advective.o advecfm2.o \
# 		ffshallw.o ffadvec.o ffeulerfm2.o \
# 		ffswfm2.o ffswfm2t.o ffeuler.o bcconv_adv.o 

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

../../doc/odocadvfm2.tex: ffadvfm2.cpp  $(dd)/advdif.tex
	$(run_odoc)

../../doc/odocburgers.tex: ffburg.cpp  $(dd)/advdif.tex
	$(run_odoc)

../../doc/odocswfm2t.tex: ffswfm2t.cpp  $(dd)/advdif.tex
	$(run_odoc)

local_dist_clean::
	-rm -f *.o

tryme: CPPFLAGS = 
tryme: secant.o tryme.o
	g++ -o tryme tryme.o secant.o 
secant.o: secant.cpp secant.h
tryme.o: tryme.cpp secant.h

secant.o: secant.cpp
	g++ -c $<

tryme.o: tryme.cpp
	g++ -c $<
