#__INSERT_LICENSE__
#$Id: Makefile,v 1.5 2001/05/26 22:18:33 mstorti Exp $

EPS = $(patsubst %.obj, %.eps, $(wildcard *.obj))
EPS_F = $(patsubst %.fig, %.eps, $(wildcard *.fig))
PDF = $(patsubst %.eps, %.pdf, $(EPS) $(EPS_F) )

.PHONY : all clean
all : $(EPS) $(EPS_F)
pdf: $(PDF)

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base

$(EPS): %.eps : %.obj
	rm -f $@
	fixobj
	objdir=`pwd` ; cd ; tgif -print -color -eps $$objdir/$<

$(EPS_F): %.eps : %.fig
	fig2dev -L eps $< $@

%.pdf: %.eps
#	ps2pdf $< $@
	epstopdf $< -o $@

local_dist_clean::
	rm -f *.eps *~ *.bak *.pdf

local_clean:: local_dist_clean
