#__INSERT_LICENSE__
#$Id: Makefile,v 1.49 2005/02/09 01:52:13 mstorti Exp $

.PHONY: all_html clean lclean all icons latex force ps readme

default: pdf

#p [in doc/Makefile]
#s Main targets

#w all the doc
all: latex index.html all_html ps pdf readme eps

PETSCFEM_DIR     = ..
include $(PETSCFEM_DIR)/Makefile.base

ODOCFILES = odocadv.tex odocadvdif.tex odocadvdif2.tex odocadvdife.tex \
	odocadve.tex odocadvfe.tex odocadvfm2.tex odocadvfs.tex \
	odocburgers.tex odocgenl.tex \
	odocns.tex odocnsfs.tex odocnsb.tex odocnse.tex odocnsw.tex odocswfm2t.tex \
	odociisd.tex odocdxhook.tex

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

dvi: petscfem.dvi

latex-petscfem: petscfem.tex
	pdflatex petscfem

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
doc_html: README.html index.html petscfem_html manual_html

README.html: readme
	cp ../README.html .

#w make html for reference only
petscfem_html: latex icons petscfem/petscfem.html 

#w make html for manual only
manual_html: 
	$(MAKE) -C manual man_html

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
	$(EPERL) $(EPERLFLAGS) $< > $@
	chmod -w $@

# Run  `$ make ltx2h_opt=-no_images' in order to not wait so much time...
petscfem/petscfem.html: odoc eps petscfem.in.tex $(TEXFILES:=.tex)
	$(MAKE) EPERLFLAGS="-d html=1" petscfem.tex
	/usr/local/latex2html/latex2html $(ltx2h_opt) \
			-init_file latex2html.init petscfem
# balance tags due to a bug in latex2html. 
	cd petscfem ; for f in *.html ; do  ../fixul.pl $$f ; done

#w Make all documentation files for program options with 'odoc.pl'
odoc: odoc_all
	$(MAKE) -C ../src odocsrc
	$(MAKE) -C ../applications/ns odocns 
	$(MAKE) -C ../applications/advective odocadv
	$(MAKE) -C ../applications/advdif odocadvdif

###--------<*>--------<*>--------<*>--------<*>--------<*>--------<*>
allcpp = $(shell find .. -name '*.cpp')
odoc_all: 
	$(ODOC) -o odocall.tex $(allcpp) &> odocall.log
	-rm odocall.log

info: 
	@$(ODOC) -e options2.texi $(allcpp) &> odoc.log
	emacs -batch -q -l make-info.el &> make-info.log
	makeinfo --no-split options3.texi
	-rm options2.texi options3.texi odoc.log make-info.log odoc.log

#w Make html home-page
index.html: index.tex
	/usr/local/latex2html/latex2html -init_file latex2html.init \
		-no_subdir -split 0 index

#w Make the README in text and html format
readme: ../README ../README.html

README.pod: README.podin ../VERSION
	-chmod +w README.pod
	$(EPERL) -P README.podin > README.pod
	-chmod -w README.pod

../README: README.pod ../VERSION
	-chmod +w $@
	cat ./banner >$@
	pod2text README.pod >> $@
	-chmod -w $@

../README.html: README.pod ../VERSION
	-chmod +w $@
	pod2html README.pod > $@
	-chmod -w $@
	-@rm pod2html-{dir,item}cache &> /dev/null

CLEAN_DIRS = manual OBJ

local_clean:: local_dist_clean
	-rm -rf *.ps.gz index petscfem icons index.html index.css		\
		version.tex *.pdf petscfem.tex					\
		contents.pl images.pl index.pl labels.pl sections.pl WARNINGS	\
		README.html options2.texi options3.texi options.info		\
		make-info.log

SWDIRS := 
local_sw:: readme 

#w Clean leaving PS and HTML
local_dist_clean::
	-putex
	-rm -rf tz*.tex .*.db \
              pod2html-{dir,item}cache petscfem.out *.bak

fm2ex.cpp: ../test/testfm2c.cpp
	$(EPERL) -B '//<' -E '//>' $< > $@

#s
