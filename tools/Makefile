#__INSERT_LICENSE__
#__INSERT_LICENSE__
.PHONY: all

scripts = odoc.pl makewhat.pl

#v $hdrmxlen += 3;
#p [in tools]
#w Builds all doc in this directory
all: $(scripts:pl=1) $(scripts:pl=txt)

PETSCFEM_DIR = ..
include $(PETSCFEM_DIR)/Makefile.base
OMKFLS = 

#w Builds man page from POD documentation in <script>.pl
#e <script>.1
%.1: %.pl
	pod2man $< > $@

#w Builds man page from POD documentation in <script>.pl
#e <script>.txt
%.txt: %.pl
	pod2text $< > $@

odoc.1: odoc.pl
odoc.txt: odoc.pl

makewhat.1: makewhat.pl
makewhat.txt: makewhat.pl

#s
