#__INSERT_LICENSE__
#$Id: Makefile,v 1.4 2001/04/01 03:10:58 mstorti Exp $

EPS = $(patsubst %.obj, %.eps, $(wildcard *.obj))
EPS_F = $(patsubst %.fig, %.eps, $(wildcard *.fig))

.PHONY : all clean
all : $(EPS) $(EPS_F)

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base

$(EPS): %.eps : %.obj
	rm -f $@
	fixobj
	objdir=`pwd` ; cd ; tgif -print -color -eps $$objdir/$<

$(EPS_F): %.eps : %.fig
	fig2dev -L ps $< $@

local_dist_clean::
	rm -f *.eps *~ *.bak

local_clean:: local_dist_clean
