#__INSERT_LICENSE__
.PHONY: all

default: run_nettest

scripts = odoc.pl makewhat.pl

#v $hdrmxlen += 3;
#p [in tools]
#w Builds all doc in this directory
all: $(scripts:pl=1) $(scripts:pl=txt) 

PETSCFEM_DIR = ..
include $(PETSCFEM_DIR)/Makefile.base
OMKFLS = 

local_sw:: 
	chmod 755 ./$(EPERL)_min ./$(EPERL)

#w all license files
licenses: license.script license.c license.tex license.cpp

#w license for makefiles and scripts
license.script: license.txt ../VERSION
	echo "#__<L>__" > $@
	echo "#__<L>__" >> $@
	perl -p -e 's/^/\#  /; s/<VERSION>/$(PETSCFEM_VERSION)/;' license.txt >> $@

#w license for C 
license.c: license.txt ../VERSION
	echo "\/.__<L>__" > $@ 
	echo "/*__<L>__*/" >> $@
	echo '/*' >>$@
	perl -p -e "s/^/  /; s/<VERSION>/$(PETSCFEM_VERSION)/;" license.txt >> $@
	echo '*/' >>$@

#w license for C++
license.cpp: license.txt ../VERSION
	echo "\/.__<L>__" > $@ 
	echo "//__<L>__" >> $@
	perl -p -e "s|^|//  |; s/<VERSION>/$(PETSCFEM_VERSION)/;" license.txt >> $@

#w license for C and C++ essence
license.tex: license.txt ../VERSION
	echo "%__<L>__" > $@
	echo "%__<L>__" >> $@
	perl -p -e 's/^/%  /; s/<VERSION>/$(PETSCFEM_VERSION)/;' license.txt >> $@

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

#w Hexasplit program splits a mesh of hexas in tetras
hexasplit.bin: hexasplit.o
	${CXX_CLINKER} -o $@ $<  $(LDFLAGS_NO_PROG_LIB)

#w 
nettest.bin: nettest.o
	${CXX_CLINKER} -o $@ $<  $(LDFLAGS_NO_PROG_LIB)

NP = $(shell procsel -n $(PETSCFEM_PROCTABLE))
# NP = $(shell procsel -n ./proctable)
#w Runs test of bandwith between all processors
run_nettest: nettest.bin
	$(MPIRUN) -np $(NP) -nolocal -machinefile ./machi.dat $<

#w Conversion form Abekas YUV format (used by DX) to standard  YUV
dxyuv: CFLAGS = -O3 -funroll-loops
dxyuv: dxyuv.cpp
	$(CXX) $(CFLAGS) $< -o $@

local_clean::
	-rm -f dxyuv

#w Generate doc for myexpect.pl
myexpect_doc: myexpect.1 myexpect.txt myexpect.html
	pod2html myexpect.pl > myexpect.html
	pod2text myexpect.pl > myexpect.txt
	pod2man myexpect.pl > myexpect.1

#s
