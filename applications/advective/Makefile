#__INSERT_LICENSE__
 # $Id: Makefile,v 1.6 2002/09/17 17:15:13 mstorti Exp $ 

.PHONY: all run lclean save compile

PROG = adv
PETSCFEM_DIR   = ../..
LIBOBJS = advective.o advecfm2.o adv.o absorb.o				\
		ffshallw.o ffadvec.o ffeulerfm2.o			\
		ffswfm2.o ffswfm2t.o bcconv_adv.o ffeuler.o bccadvfm2.o	\
		compdate.o
PETSCFEM_MODULE = Advective
MYOBJS = 
TAGFLAGS = $(TAGINCL)="../../src/TAGS"
#TAGFLAGS = --include="../../src/TAGS"

SRCS := $(MYOBJS:.o=.cpp)

include $(PETSCFEM_DIR)/Makefile.base

compile: $(PROG).bin

%.cppi: %.cpp
	g++ $(CPPFLAGS) -I /usr/local/mpich-1.2.0/include/ \
		-I/home/mstorti/PETSC/petsc-2.0.24/bmake/linux/ \
		-I/home/mstorti/PETSC/petsc-2.0.24/include/ -E $< > $@

ffeuler.o: ffeuler.cpp
	-${CC} -c ${COPTFLAGS} ${CFLAGS} ${CCPPFLAGS} -O0 $<

# documentation of options with the script `odoc.pl'
odocadv: ../../doc/odocadv.tex  ../../doc/odocadve.tex \
	../../doc/odocadvfe.tex ../../doc/odocadvfs.tex 

../../doc/odocadv.tex: adv.cpp
	$(run_odoc)

../../doc/odocadve.tex: advecfm2.cpp
	$(run_odoc)

../../doc/odocadvfe.tex: ffeulerfm2.cpp
	$(run_odoc)

../../doc/odocadvfs.tex: ffswfm2.cpp
	$(run_odoc)

