#__INSERT_LICENSE__
#$Id: Makefile,v 1.21 2001/08/06 01:06:10 mstorti Exp $

.PHONY: html clean lclean all icons latex force ps readme

#p [in doc/Makefile]
#s Main targets

#w all the doc
all: latex index.html html ps pdf readme eps

PETSCFEM_DIR     = ..
include $(PETSCFEM_DIR)/Makefile.base

ODOCFILES = odocadv.tex odocadvdif.tex odocadvdife.tex \
	odocadve.tex odocadvfe.tex odocadvfm2.tex odocadvfs.tex \
	odocburgers.tex odocgenl.tex \
	odocns.tex odocnsb.tex odocnse.tex odocnsw.tex odocswfm2t.tex \
	odociisd.tex

TEXFILES = prev1 layout thash idmap compprof advec \
		ns tests fastmat2 $(ODOCFILES:.tex=)

uncomment:
	perl -pi.bak -e 's/^%+\s*(.*%__UNCOMMENT_IN_FINAL_VERSION__)$$/$$1/;' \
			-e 's/^([^%].*%__COMMENT_IN_FINAL_VERSION__)$$/%$$1/;' \
			petscfem.in.tex

#w makes petscfem.ps.gz, manual.ps.gz
ps: petscfem.ps.gz manual.ps.gz

#w make PDF
pdf: petscfem.pdf

#w Makes PS version of the guide
#e petscfem.ps.gz

#w Makes PS version of the reference
#e manual.ps.gz

.INTERMEDIATE: manual.ps petscfem.ps
petscfem.ps: petscfem.dvi
	dvips -o $@ petscfem

petscfem.dvi: *.tex
	$(MAKE) latex

petscfem.pdf: latex
	$(MAKE) -C OBJ pdf
	pdflatex petscfem

manual.ps.gz: force
	$(MAKE) -C manual ps

icons:
	if [ ! -d icons ] ; then uudecode < icons.tgz.uue ; \
		tar -xzvf icons.tgz ; rm icons.tgz ; fi

#w make html for reference and manual
html: index.html petscfem_html manual_html

#w make html for reference only
petscfem_html: latex icons petscfem/petscfem.html 

#w make html for manual only
manual_html: 
	$(MAKE) -C manual html

#s Other options

#w process LaTeX on reference
latex: force eps odoc
	$(MAKE) EPERLFLAGS="-d html=0" petscfem.tex
	for j in 1 2 3 4 ; \
	do echo "On run $$j -------------------------------" ; \
		latex petscfem ; \
		bibtex petscfem ; \
		makeindex petscfem.idx ; \
		if ! grep 'LaTeX Warning: .* Rerun ' petscfem.log ; then break ; fi ; \
	done

#w make .eps version of the figures
eps:
	$(MAKE) -C OBJ all
	if [ ! -d petscfem/OBJ ] ; \
		then mkdir petscfem ; cd petscfem ; ln -s ../OBJ . ; \
	fi

petscfem.tex: petscfem.in.tex force
	if [ -f $@ ] ; then rm -f $@ ; fi
	eperl $(EPERLFLAGS) $< > $@
	chmod -w $@

# Run  `$ make ltx2h_opt=-no_images' in order to not wait so much time...
petscfem/petscfem.html: odoc eps petscfem.in.tex $(TEXFILES:=.tex)
	$(MAKE) EPERLFLAGS="-d html=1" petscfem.tex
#	latex2html $(ltx2h_opt) -prefix pfem -init_file latex2html.init 
	latex2html $(ltx2h_opt) -init_file latex2html.init \
		-split 2 -reuse 2 petscfem

#w Make all documentation files for program options with 'odoc.pl'
odoc:
	$(MAKE) -C ../src odocsrc
	$(MAKE) -C ../applications/ns odocns
	$(MAKE) -C ../applications/advective odocadv
	$(MAKE) -C ../applications/advdif odocadvdif

#w Make html home-page
index.html: index.tex
	latex2html -init_file latex2html.init -no_subdir -split 0 index

#w Make the README in text and html format
readme: ../README ../README.html

README.pod: README.podin
# 	-chmod +w README.pod
# 	fversion -is $(TARFILE).tara,v README.podin README.pod
# 	-chmod -w README.pod
	-chmod +w README.pod
	eperl -P README.podin > README.pod
	-chmod -w README.pod

../README: README.pod
	-chmod +w $@
	cat ./banner >$@
	pod2text $< >> $@
	-chmod -w $@

../README.html: README.pod
	-chmod +w $@
	pod2html $< > $@
	-chmod -w $@
	-@rm pod2html-{dir,item}cache &> /dev/null

CLEAN_DIRS = manual OBJ

local_clean:: local_dist_clean
	-rm -rf *.ps.gz index petscfem icons index.html index.css \
		version.tex *.pdf petscfem.tex \
		contents.pl images.pl index.pl labels.pl sections.pl 

#w Clean leaving PS and HTML
local_dist_clean::
	-putex
	-rm -rf tz*.tex .*.db \
              pod2html-{dir,item}cache petscfem.out *.bak

fm2ex.cpp: ../test/testfm2c.cpp
	eperl -B '//<' -E '//>' $< > $@

#s
