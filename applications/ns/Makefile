#__INSERT_LICENSE__
 # $Id: Makefile,v 1.23 2002/02/13 19:42:08 mstorti Exp $ 

.PHONY: all run lclean save compile

PROG = ns
PETSCFEM_DIR   = ../..
#MYOBJS = ns.o fracstep.o nsi_tet.o

LIBOBJS = ns.o nsitetlesfm2.o					\
		bccnsfm2.o wall.o nsvol.o walldata.o		\
		nsilesther.o bccnstherfm2.o nsikeps.o		\
		wallke.o nonlr.o lagmul.o nssup.o nssupr.o	\
		nsid.o nsres.o adaptor.o elast.o adappg.o	\
		qharm.o

MYOBJS = 
#MYOBJS = ns.$(opt_level).o 

%.g.o: %.o
#	If object does not have the appropriate debugging level recompile
#	[ OK. This is very tricky... I don't know the portability of this ]
	@if objdump -g $< 2>&1 | grep 'no recognized debugging information'	\
		> /dev/null ; then						\
		obj_optl='O' ;							\
		else								\
		obj_optl='g' ;							\
	fi ;									\
	if [ ! "xg"="x$$obj_optl" ] ; then					\
		echo "OptLevelCheck: object $< does not match optimization "	\
				"level BOPT: $(BOPT)" ;				\
		rm $< ;								\
		$(MAKE) $< ;							\
	fi ;									\
	mv -fv $< $@

%.O.o: %.o
#	If object does not have the appropriate debugging level recompile
#	[ OK. This is very tricky... I don't know the portability of this ]
	@if objdump -g $< 2>&1 | grep 'no recognized debugging information'	\
		> /dev/null ; then						\
		obj_optl='O' ;							\
		else								\
		obj_optl='g' ;							\
	fi ;									\
	if [ ! "xO"="x$$obj_optl" ] ; then					\
		echo "OptLevelCheck: object $< does not match optimization "	\
				"level BOPT: $(BOPT)" ;				\
		rm $< ;								\
		$(MAKE) $< ;							\
	fi ;									\
	mv -fv $< $@

TAGFLAGS = $(TAGINCL)="../../src/TAGS"
#TAGFLAGS = --include="../../src/TAGS"

# for makedepend
SRCS := $(MYOBJS:.$(opt_level).o=.cpp)

include $(PETSCFEM_DIR)/Makefile.base

PROG_LIB = lib$(PROG)_$(opt_level).a

$(PROG_LIB): $(PROG_LIB)($(LIBOBJS))
	ranlib $@

$(PROG_LIB)(%.o): %.o
	ar ru $@ $<

$(PROG) compile: $(PROG).bin

%.cppi: %.cpp
	g++ $(CPPFLAGS) -I /usr/local/mpich-1.2.0/include/		\
		-I/home/mstorti/PETSC/petsc-2.0.24/bmake/linux/		\
		-I/home/mstorti/PETSC/petsc-2.0.24/include/ -E $< > $@

# documentation of options with the script `odoc.pl'
odocns: ../../doc/odocns.tex ../../doc/odocnse.tex \
	../../doc/odocnsb.tex  ../../doc/odocnsw.tex 

../../doc/odocns.tex: ns.cpp 
	$(run_odoc)

../../doc/odocnse.tex: nsitetlesfm2.cpp 
	$(run_odoc)

../../doc/odocnsb.tex: bccnsfm2.cpp 
	$(run_odoc)

../../doc/odocnsw.tex: wall.cpp 
	$(run_odoc)

noafs:
	@echo "making symbolic links to NO Algebraic FS formulation... "
	@ls -l ns.cpp nsitetlesfm2.cpp ; \
	echo -n "re-link to .noafs files? (y/n) > " ; \
	read answer ; \
	if [ -L ns.cpp -a -L nsitetlesfm2.cpp -a "$$answer" = "y" ] ; \
		then ln -sf ns.cpp.noafs ns.cpp ; \
		ln -sf nsitetlesfm2.cpp.noafs nsitetlesfm2.cpp ; \
	fi ; 
	-@rm ns.o nsitetlesfm2.o

afs:
	@echo "making symbolic links to Algebraic FS formulation... "
	@ls -l ns.cpp nsitetlesfm2.cpp ; \
	echo -n "re-link to .afs files? (y/n) > " ; \
	read answer ; \
	if [ -L ns.cpp -a -L nsitetlesfm2.cpp -a "$$answer" = "y" ] ; \
		then ln -sf ns.cpp.afs ns.cpp ; \
		ln -sf nsitetlesfm2.cpp.afs nsitetlesfm2.cpp ; \
	fi ; 
	-@rm ns.o nsitetlesfm2.o

#lclean: general_clean
#	cat /dev/null >makefile.d
#	-rm *.o $(OBJS) $(PROG).bin octave-core *.cppi outvector.sal makefile.d
#	$(MAKE) clean

# File containing dependencies computed by makedepend
#
