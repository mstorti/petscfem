# mode: -*- makefile -*-
#__INSERT_LICENSE__
#$Id: Makefile,v 1.48 2003/03/17 19:02:49 mstorti Exp $

SHELL = /bin/bash

.PHONY: all run lclean save libpetscfem ns adv laplace doc newdepend tags	\
		sw startwork fm2new sync_version applications			\
		ns_O ns_g

APPS = adv advdif ns laplace
APPDIRS = advective advdif ns laplace
APPDIRS := $(patsubst %,applications/%,$(APPDIRS))

CLEAN_DIRS := src $(APPDIRS) doc test tools dx

SRCDIRS := src $(APPDIRS) test 

SRCS := 

DEPEND_DIRS := $(SRCDIRS)

SWDIRS := test tools

#p [in Makefile]

#s Main targets
#w Makes doc, library and applications. No cleaning 
all: sw doc pflib $(APPS)

#w Builds all necessary things after checking out a version
#w from the CVS repository
local_sw:: sync_version depend tags 
	cd tools; ln -sf hexenco.pl ident2iso ; ln -sf hexenco.pl iso2ident
	$(MAKE) -C doc readme
# make scripts executable
	chmod 755 ./make/appatch ./make/mkpatch ./make/mkvers \
		./src/insdeb.pl ./test/runtests.pl ./test/turbchan/verify.pl \
		./tools/coall ./tools/checktag \
		./tools/insert_license.pl ./tools/makeltag \
		./tools/maketag ./tools/makewhat.pl ./tools/myexpect.pl \
		./tools/odoc.pl ./tools/petscload.pl ./tools/pfcpp \
		./doc/manual/vrfdocpp.pl ./src/insdeb.pl ./doc/fixul.pl

local_clean::
	cd tools ; rm -f ident2iso ; rm -f iso2ident

#w Builds existng applications
applications: 
	$(MAKE) $(APPS)
	for dir in $(APPDIRS) ; do $(MAKE) -C $$dir distclean ; done

#w Builds a package + doc + applications 
bin_distrib: license sw finaldoc pflib applications

#w Builds a package + doc 
distrib: license sw finaldoc pflib applications

C_LICENSE_FILE = tools/license.c
SCRIPT_LICENSE_FILE = tools/license.script
TEX_LICENSE_FILE = tools/license.tex

#w Includes license in all .cpp and .h files
license:
	$(MAKE) -C tools licenses
	find . -type f -name '*.cpp' -print0 | xargs -0 -n 1 -e \
		$(INSERT_LICENSE) $(C_LICENSE_FILE)
	find . -type f -name '*.h' -print0 | xargs -0 -n 1 -e \
		$(INSERT_LICENSE) $(C_LICENSE_FILE)
	find . -type f -name '*.y' -print0 | xargs -0 -n 1 -e \
		$(INSERT_LICENSE) $(C_LICENSE_FILE)
	find . -type f -name 'Makefile*' -print0 | xargs -0 -n 1 -e \
		$(INSERT_LICENSE) $(SCRIPT_LICENSE_FILE)
	find . -type f -name '*.m' -print0 | xargs -0 -n 1 -e \
		$(INSERT_LICENSE) $(SCRIPT_LICENSE_FILE)
	find . -type f -name '*.pl' -print0 | xargs -0 -n 1 -e \
		$(INSERT_LICENSE) $(SCRIPT_LICENSE_FILE)
	find . -type f -name '*.tex' -print0 | xargs -0 -n 1 -e \
		$(INSERT_LICENSE) $(TEX_LICENSE_FILE)

unlicense:
	find . -type f -name '*.cpp' -print0 | xargs -0 -n 1 -e \
		$(INSERT_LICENSE) -u $(C_LICENSE_FILE)
	find . -type f -name '*.h' -print0 | xargs -0 -n 1 -e \
		$(INSERT_LICENSE) -u $(C_LICENSE_FILE)
	find . -type f -name '*.y' -print0 | xargs -0 -n 1 -e \
		$(INSERT_LICENSE) -u $(C_LICENSE_FILE)
	find . -type f -name 'Makefile*' -print0 | xargs -0 -n 1 -e \
		$(INSERT_LICENSE) -u $(SCRIPT_LICENSE_FILE)
	find . -type f -name '*.m' -print0 | xargs -0 -n 1 -e \
		$(INSERT_LICENSE) -u $(SCRIPT_LICENSE_FILE)
	find . -type f -name '*.pl' -print0 | xargs -0 -n 1 -e \
		$(INSERT_LICENSE) -u $(SCRIPT_LICENSE_FILE)
	find . -type f -name '*.tex' -print0 | xargs -0 -n 1 -e \
		$(INSERT_LICENSE) -u $(TEX_LICENSE_FILE)

#w Uncomments \input lines for partial processing and then makes doc
finaldoc:
	$(MAKE) -C doc uncomment
	$(MAKE) doc

#w Builds the doc
doc:
	$(MAKE) sync_version
	$(MAKE) -C doc all distclean

#w Builds the library
pflib:
	$(MAKE) -C src compile distclean

TAGDIRS = src $(APPDIRS)
#s other targets
#w Builds/refresh Emacs tags tables
tags: 
	for dir in $(TAGDIRS) ; do $(MAKE) -C $$dir TAGS ; done

PETSCFEM_DIR     = .
include $(PETSCFEM_DIR)/Makefile.base
LOCDIR           = $(PWD)

DIRS = doc manual src ns advective tryme laplace 

#w Makes tests
tests:
	$(MAKE) -C test tests

#----<*>----<*>----<*>----<*>----<*>----<*>----<*>----<*>----
# APPLICATIONS


#w Builds the Navier Stokes module
ns: libpetscfem 
	$(MAKE) -C applications/ns ns$(osfx).bin

#w Make NS module (debugger and optimized versions)
ns_all: ns_g ns_O

#w Make NS module (debugger version)
ns_g: 
	$(MAKE) BOPT=g_c++ ns

#w Make NS module (optimized version)
ns_O:
	$(MAKE) BOPT=O_c++ ns

#w Make Advdif module (O_c++ and g_c++)
advdif_all: 
	$(MAKE) BOPT=g_c++ advdif
	$(MAKE) BOPT=O_c++ advdif

#w Builds th advective systems module (Euler, shallow water)
adv: libpetscfem
	$(MAKE) -C applications/advective adv$(osfx).bin

#w Builds the Advective/diffusive systems module (NS-compresible,
#w 		(shallow water+diffusive and turbulent terms, 
#w               linear advection diffusion, burgers
advdif: libpetscfem
	$(MAKE) -C applications/advdif advdif$(osfx).bin

#w Builds the Laplace module
laplace: libpetscfem
	$(MAKE) -C applications/laplace laplace$(osfx).bin

#----<*>----<*>----<*>----<*>----<*>----<*>----<*>----<*>---- 
%.cppi: %.cpp
	g++ -E $(CCPPFLAGS) $< > $@ ; chmod u-w $@

#w Resyncs some administrative files with the current version number.
sync_version: 	
	@version=`cat VERSION` ;					\
	echo "Creating doc/version.tex" ;				\
	echo "\\def\\petscfemversion{$$version}" > doc/version.tex ;	\
	echo "Creating src/version.cpp" ;				\
	$(MAKE) -C src version.cpp
	$(MAKE) -C doc readme

#w This is somewhat obsolet. Now uses CVS
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

#w Count lines in `.cpp' and `.h' files
line_count:
	@echo -n 'Lines, word, chars in source files: ' ; \
	cat `find . -name '*.cpp'` `find . -name '*.h'` \
	`find . -name '*.c'` `find . -name '*.f'` | wc 

#w Makes a new release
tag:
	./tools/maketag

#w Makes a new (light) release
ltag:
	./tools/makeltag save.log

torture:
	$(MAKE) -C test torture

#s
