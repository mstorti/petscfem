# $Id: Makefile,v 1.4 2001/01/04 20:06:02 mstorti Exp $ 

.PHONY: all run lclean save libpetscfem ns adv laplace doc newdepend tags \
		sw startwork fm2new

all: newdepend tags adv ns laplace 

#sw: newdepend tags ns
sw: newdepend tags advdif

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

CLEANDIRS = src applications/laplace applications/advective \
			applications/ns applications/advdif \
			doc test run run/LES run/algebfs

SRCDIRS = src applications/ns applications/advective applications/advdif \
		applications/laplace test 

DEPEND_DIRS = $(SRCDIRS)

depend:
#	for dir in $(SRCDIRS) ; do touch $$dir/makefile.d ; done ; \
#		for dir in $(SRCDIRS) ; do $(MAKE) -C $$dir depend ; done
	for dir in $(SRCDIRS) ; do $(MAKE) -C $$dir depend ; done 

sync_version: 	
	$(MAKE) -C src version.cpp
	$(MAKE) -C doc readme

save:
	$(MAKE) sync_version
	if [ -f $(TARFILE).tara,v ]  ; \
		then co  -ko -l $(TARFILE).tara ; \
		mv -f $(TARFILE).tara $(TARFILE).tara.old ; \
		fi 
	tar -cvf $(TARFILE).tar -X make/exclude -C .. `thisdir.pl`
	tar2tara < $(TARFILE).tar > $(TARFILE).tara ; rm $(TARFILE).tar
	@echo "Old and new tara files: "
	@ls -l $(TARFILE).tara.old $(TARFILE).tara
	@echo "Differences (lines/words/chars):"
	@diff $(TARFILE).tara.old $(TARFILE).tara | wc 
	@echo -n "Proceed to commit? (y/n) " 
	@read answer ; \
	if [ $$answer = "y" ] ; then \
		echo "Checking in $(TARFILE) ..." ; \
		ci $(TARFILE).tara ; \
		echo "Generating history ..." ; \
		rlog $(TARFILE).tara > petscfem_hist.txt ; \
		pwd ; \
		echo "Getting version ..." ; \
		$(MAKE) sync_version ; \
	else rm  $(TARFILE).tara ; \
	fi ;\
	if [ -e $(TARFILE).tara.old ] ; then rm $(TARFILE).tara.old ; fi



