# doc/manual/Makefile for DOC++
# Copyright (c) 1999 Dragos Acostachioaie

.PHONY: html tex

#DOCFILES = ../src/{idmap,fem,utils}.h
DOCF = doc.dxx
DOCXX = doc++

all:	html

html:
	$(DOCXX) -p -d html $(DOCF)

tex:

ps: ../manual.ps

../manual.ps: manual.dvi
	dvips -o ../manual.ps manual.dvi

manual.dvi: manual.tex
	-latex \\nonstopmode\\input manual
	-latex \\nonstopmode\\input manual

manual.tex:
	$(DOCXX) -d tex -v -t -p -o manual.tex -ep a4wide $(DOCF)
# Fix some error in DOC++
	perl -p -i.bak \
		-e "s/classTrandom_popsetT/class T random\\\\_pop(set\\\\<T\\\\> \\\\&)/" \
				manual.tex

clean:
	$(MAKE) -C html -f ../Makefile html_clean 
	-rm manual.tex manual.tex.bak
	-putex

html_clean:
	rm -rf *.html *.tex *.ps *.dvi *.gif *.class \
		gifs.db *.bak *.aux *.eps *.pnm *~ *.log
