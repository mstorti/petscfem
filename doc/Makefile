#$Id: Makefile,v 1.6 2001/01/09 16:58:51 mstorti Exp $

.PHONY: html clean lclean all icons latex force ps readme

all: latex index.html html ps readme eps

PETSCFEM_DIR     = ..
include $(PETSCFEM_DIR)/Makefile.base

ODOCFILES = odocadv.tex odocadve.tex odocadvfe.tex odocadvfs.tex \
	odocnsw.tex odocns.tex odocnsb.tex odocnse.tex

TEXFILES = prev1 layout thash idmap compprof advec \
		ns tests fastmat2 $(ODOCFILES:.tex=)

ps: petscfem.ps manual.ps

petscfem.ps: petscfem.dvi
	dvips -o $@ petscfem

petscfem.dvi: *.tex
	$(MAKE) latex

manual.ps: force
	$(MAKE) -C manual ps

icons:
	if [ ! -d icons ] ; then uudecode < icons.tgz.uue ; \
		tar -xzvf icons.tgz ; rm icons.tgz ; fi

html: index.html petscfem_html manual_html

petscfem_html: latex icons petscfem/petscfem.html 

manual_html: 
	$(MAKE) -C manual html

latex: force eps odoc
	latex petscfem
	-bibtex petscfem
	-makeindex petscfem
	latex petscfem
	latex petscfem

eps:
	$(MAKE) -C OBJ all
	if [ ! -d petscfem/OBJ ] ; \
		then mkdir petscfem ; cd petscfem ; ln -s ../OBJ . ; \
	fi


petscfem/petscfem.html: odoc eps petscfem.tex petscfem*.tex $(TEXFILES:=.tex)
	latex2html $(ltx2h_opt) -prefix pfem -init_file latex2html.init \
		-split 2 -reuse 2 petscfem

odoc:
	$(MAKE) -C ../applications/ns odocns
	$(MAKE) -C ../applications/advective odocadv
	$(MAKE) -C ../applications/advdif odocadvdif

#odocns.tex odocnse.tex odocnsb.tex odocnsw.tex: 
#odocadv.tex odocadve.tex odocadvfe.tex odocadvfs.tex: 

index.html: index.tex
	latex2html -init_file latex2html.init -no_subdir -split 0 index

readme: ../README ../README.html

README.pod: README.podin $(TARFILE).tara,v
# 	-chmod +w README.pod
# 	fversion -is $(TARFILE).tara,v README.podin README.pod
# 	-chmod -w README.pod
	-chmod +w README.pod
	eperl -P README.podin > README.pod
	-chmod -w README.pod

../README: README.pod
	-chmod +w $@
	pod2text $< > $@
	-chmod -w $@

../README.html: README.pod
	-chmod +w $@
	pod2html $< > $@
	-chmod -w $@
	-@rm pod2html-{dir,item}cache &> /dev/null

local_clean::
	$(MAKE) -C manual clean
	$(MAKE) -C OBJ clean
	-putex
	-rm tz*.tex index.html index.css .*.db
	-rm -rf index
	-rm -rf petscfem
	-rm -rf icons pod2html-{dir,item}cache

fm2ex.cpp: ../test/testfm2c.cpp
	eperl -B '//<' -E '//>' $< > $@
