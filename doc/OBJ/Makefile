#$Id: Makefile,v 1.1 2000/12/28 12:54:42 mstorti Exp $

EPS = $(patsubst %.obj, %.eps, $(wildcard *.obj))
EPS_F = $(patsubst %.fig, %.eps, $(wildcard *.fig))

.PHONY : all clean
all : $(EPS) $(EPS_F)

$(EPS): %.eps : %.obj
	rm -f $@
	fixobj
	objdir=`pwd` ; cd ; tgif -print -color -eps $$objdir/$<

$(EPS_F): %.eps : %.fig
	fig2dev -L ps $< $@
clean: 
	rm -f *.eps *~ *.bak
