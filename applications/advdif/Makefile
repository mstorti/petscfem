 # $Id: Makefile,v 1.8 2001/01/23 13:25:16 mstorti Exp $ 

.PHONY: all run lclean save compile

PROG = advdif
PETSCFEM_DIR   = ../..
MYOBJS = advdife.o advdif.o ffadvfm2.o  bccadvfm2.o  ffburg.o
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

local_dist_clean::
	-rm -f *.o
