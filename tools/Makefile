#  [-begin-license-]
#  This  file  belongs  to  the  PETSc  -  FEM  package,  a  library  and
#  application  suite oriented  to  the Finite  Element  Method based  on
#  PETSc.   Copyright  (C)  1999-2001,  Mario  Alberto  Storti,  Norberto
#  Marcelo Nigro, Centro Internacional de Metodos Numericos en Ingenieria
#  (CIMEC-Argentina),  Universidad Nacional del  Litoral (UNL-Argentina),
#  Consejo   Nacional   de   Investigaciones   Cientificas   y   Tecnicas
#  (UNL-Argentina).
#  
#  This program is  free software; you can redistribute  it and/or modify
#  it under the  terms of the GNU General Public  License as published by
#  the Free Software Foundation; either  version 2 of the License, or (at
#  your option) any later version.
#  
#  This program  is distributed in the  hope that it will  be useful, but
#  WITHOUT   ANY  WARRANTY;   without  even   the  implied   warranty  of
#  MERCHANTABILITY  or FITNESS  FOR A  PARTICULAR PURPOSE.   See  the GNU
#  General Public License for more details.
#  
#  You  should have received  a copy  of the  GNU General  Public License
#  along  with  this  program;  if   not,  write  to  the  Free  Software
#  Foundation, Inc.,  59 Temple Place, Suite 330,  Boston, MA 02111-1307,
#  USA.
#  [-end-license-]


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
	echo $(VERSION)
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
