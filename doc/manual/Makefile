#__INSERT_LICENSE__
# doc/manual/Makefile for DOC++
# Copyright (c) 1999 Dragos Acostachioaie

.PHONY: html ps

#DOCFILES = ../src/{idmap,fem,utils}.h
DOCF = doc.dxx
DOCXX = doc++

all: html ps

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base

#w generate HTML
html:
	$(DOCXX) -p -a -ep mstorti -d html $(DOCF)

#w Generate PS
ps: ../manual.ps.gz

.INTERMEDIATE: ../manual.ps 
../manual.ps.gz: ../manual.ps
.SECONDARY: manual.tex

../manual.ps: manual.dvi
	dvips -o $@ $<

#w TeX compile manual.tex
manual.dvi: manual.tex
	-latex \\nonstopmode\\input manual
	-latex \\nonstopmode\\input manual

#w make the tex from the doc++ comments
manual.tex:
	$(DOCXX) -d tex -v -t -p -o manual.tex -ep a4wide $(DOCF)
# Fix some error in DOC++
	perl -p -i.bak \
		-e "s/classTrandom_popsetT/class T random\\\\_pop(set\\\\<T\\\\> \\\\&)/" \
				manual.tex

local_dist_clean::
	-putex
	-rm -f manual.tex.bak

local_clean:: local_dist_clean
	-rm -f manual.tex manual.pdf
	if [ -d html ] ; then rm -rf html ; fi

html_clean:
	rm -rf *.html *.tex *.ps *.dvi *.gif *.class \
		gifs.db *.bak *.aux *.eps *.pnm *~ *.log
