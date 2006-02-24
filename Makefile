# mode: -*- makefile -*-
#__INSERT_LICENSE__
#$Id: Makefile,v 1.63.4.1 2006/02/24 12:26:57 rodrigop Exp $

SHELL = /bin/bash

.PHONY: all run lclean save libpetscfem ns adv laplace doc newdepend tags	\
		sw startwork fm2new sync_version applications apps 		\
		ns_O ns_g

APPS = ns advdif adv laplace
APPDIRS = ns advdif advective laplace
APPDIRS := $(patsubst %,applications/%,$(APPDIRS))

CLEAN_DIRS := src $(APPDIRS) doc test tools dx

SRCDIRS := src $(APPDIRS) test 

SRCS := 

DEPEND_DIRS := $(SRCDIRS)

SWDIRS := tools test doc

#p [in Makefile]

#s Main targets
#w Makes doc, library and applications. No cleaning 
all: sw doc pflib $(APPS)

#w Builds all necessary things after checking out a version
#w from the CVS repository
local_sw:: 
# 	First of all make scripts executable
	chmod 755 tools/eperl tools/eperl_min ./make/appatch 		\
		./make/mkpatch ./make/mkvers				\
		./src/insdeb.pl ./test/runtests.pl 			\
		./test/turbchan/verify.pl				\
		./tools/coall ./tools/checktag ./tools/eperl		\
		./tools/insert_license.pl ./tools/makeltag		\
		./tools/maketag ./tools/makewhat.pl ./tools/myexpect.pl	\
		./tools/odoc.pl ./tools/petscload.pl ./tools/pfcpp	\
		./doc/manual/vrfdocpp.pl ./src/insdeb.pl ./doc/fixul.pl
	$(MAKE) sync_version depend 
	cd tools; ln -sf hexenco.pl ident2iso ; ln -sf hexenco.pl iso2ident

local_clean::
	cd tools ; rm -f ident2iso ; rm -f iso2ident

#w Builds existing applications
apps applications: 
	$(MAKE) $(APPS)
	for dir in $(APPDIRS) ; do $(MAKE) -C $$dir distclean ; done

#w Builds existing applications with both opt. levels (g_c++ and O_c++)
apps_all: 
	$(MAKE) BOPT=g_c++ apps
	$(MAKE) BOPT=O_c++ apps

#w Make all `new' (NS and Advdif) applications
napps_all: 
	$(MAKE) advdif_all ns_all

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
	@echo -n 'Lines in source files: ' ;	\
	cat					\
	`find . -name '*.cpp'`			\
	`find . -name '*.h'`			\
	`find . -name '*.c'`			\
	`find . -name '*.f'`			\
	`find . -name '*.tex'`			\
	`find . -name '*.pl'`			\
	`find . -name '*.m'`			\
	| wc -l

#w Makes a new release
tag:
	./tools/maketag

#w Makes a new (light) release
ltag:
	./tools/makeltag save.log

torture:
	$(MAKE) -C test torture

#w Cleans all optimized libraries and binaries
clean_O:
	-rm -f `find . -name '*_O.a'` `find . -name '*_O.bin'` 

clean_g:
	-rm -f `find . -name '*_g.a'` `find . -name '*_g.bin'` 

#w Updates working directory
sync:
	cvs up .

#s
pp:
	echo $(BS_INCLUDE)
