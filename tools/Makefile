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

#w all license files
licenses: license.script license.c license.tex

#w license for makefiles and scripts
license.script: license.txt
	echo "#__<L>__" > $@
	echo "#__<L>__" >> $@
	perl -p -e 's/^/\#  /; s/<VERSION>/$(PETSCFEM_VERSION)/;' license.txt >> $@
	echo >>$@
	echo >>$@
#	echo "#--------------------------------" >>$@

#w license for C and C++
license.c: license.txt
	echo "\/.__<L>__" > $@ 
	echo "/*__<L>__*/" >> $@
	echo '/*' >>$@
	perl -p -e "s/^/  /; s/<VERSION>/$(PETSCFEM_VERSION)/;" license.txt >> $@
	echo '*/' >>$@
	echo >>$@
	echo >>$@
#	echo "/*------------------------------*/" >>$@

#w license for C and C++
license.tex: license.txt
	echo "%__<L>__" > $@
	echo "%__<L>__" >> $@
	perl -p -e 's/^/%  /; s/<VERSION>/$(PETSCFEM_VERSION)/;' license.txt >> $@
	echo >>$@
	echo >>$@
#	echo "%--------------------------------" >>$@

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
