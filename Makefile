# $Id: Makefile,v 1.12 2001/01/10 22:40:01 mstorti Exp $ 
SHELL = /bin/bash

.PHONY: all run lclean save libpetscfem ns adv laplace doc newdepend tags \
		sw startwork fm2new sync_version

all: sync_version depend tags adv advdif ns laplace 

#sw: newdepend tags ns
sw: sync_version depend tags

fm2new:
	cd src ; ln -sf fmat2ep.cpp.new fmat2ep.cpp
	$(MAKE) fm2

fm2bck:
	cd src ; ln -sf fmat2ep.cpp.bck fmat2ep.cpp
	$(MAKE) fm2

fm2:
	$(MAKE) -C src -W fmat2ep.cpp compile
	$(MAKE) ns

TAGDIRS = src applications/advective applications/ns applications/laplace \
		applications/advdif
tags: 
	for dir in $(TAGDIRS) ; do $(MAKE) -C $$dir TAGS ; done

PETSCFEM_DIR     = .
include $(PETSCFEM_DIR)/Makefile.base
LOCDIR           = $(PWD)

DIRS = doc manual src ns advective tryme laplace 

# the library
libpetscfem:
	$(MAKE) -C src compile

# applications

# Navier Stokes
ns: libpetscfem
	$(MAKE) -C applications/ns compile

# Advective systems Euler/Shallow water
adv: libpetscfem
	$(MAKE) -C applications/advective compile

# Advective/diffusive systems (NS-compresible)/
# 		(Shallow water+diffusive and turbulent terms)
advdif: libpetscfem
	$(MAKE) -C applications/advdif compile

# Laplace equation
laplace: libpetscfem
	$(MAKE) -C applications/laplace compile

%.cppi: %.cpp
	g++ -E $(CCPPFLAGS) $< > $@ ; chmod u-w $@

CLEAN_DIRS = src applications/laplace applications/advective \
			applications/ns applications/advdif \
			doc test run run/LES run/algebfs run/sqcav

SRCDIRS = src applications/ns applications/advective applications/advdif \
		applications/laplace test 

SRCS = 

DEPEND_DIRS = $(SRCDIRS)

sync_version: 	
	@version=`cat VERSION` ;					\
	echo "Creating doc/version.tex" ;				\
	echo "\\def\\petscfemversion{$$version}" > doc/version.tex

save:
	$(MAKE) sync_version
	if [ -f $(TARFILE).tara,v ]  ;				\
		then co  -ko -l $(TARFILE).tara ;		\
		mv -f $(TARFILE).tara $(TARFILE).tara.old ;	\
		fi 
	tar -cvf $(TARFILE).tar -X make/exclude -C .. `thisdir.pl`
	tar2tara < $(TARFILE).tar > $(TARFILE).tara ; rm $(TARFILE).tar
	@echo "Old and new tara files: "
	@ls -l $(TARFILE).tara.old $(TARFILE).tara
	@echo "Differences (lines/words/chars):"
	@diff $(TARFILE).tara.old $(TARFILE).tara | wc 
	@echo -n "Proceed to commit? (y/n) " 
	@read answer ;								\
	if [ $$answer = "y" ] ; then						\
		echo "Checking in $(TARFILE) ..." ;				\
		ci $(TARFILE).tara ;						\
		echo "Generating history ..." ;					\
		rlog $(TARFILE).tara > petscfem_hist.txt ;			\
		pwd ;								\
		echo "Getting version ..." ;					\
		$(MAKE) sync_version ;						\
	else rm  $(TARFILE).tara ;						\
	fi ;									\
	if [ -e $(TARFILE).tara.old ] ; then rm $(TARFILE).tara.old ; fi

tag:
	echo shell: $(SHELL)
	@echo Verify that current directory is OK...
	echo n | cvs release .
	@echo -n "Continue? (y/n) > " ;		\
	read answer ;				\
	if [ ! $$answer = "y" ] ; then		\
		exit ;				\
	fi
	@echo Last tags:
	@grep "^tag: " save.log | tail
	@echo -n "Enter new tag: > " ;					\
	read newtag ;							\
	newtag_=`echo $$newtag | perl -pe 's/\-/--/g; s/\./-/g;'` ;	\
	echo "encoded tag: $$newtag_" ;					\
	echo "tag: $$newtag on `date`, by "				\
		"`whoami` in `hostname -f`" >> save.log ;		\
	echo $$newtag > VERSION ;					\
	echo "Proceed to tag files (y/n) > " ;				\
	read answer ;							\
	if [ $$answer = "y" ] ; then					\
		cvs tag $$newtag_ . ;					\
	fi
