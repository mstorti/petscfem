 # $Id: Makefile,v 1.1.1.1 2000/12/28 12:54:43 mstorti Exp $ 

.PHONY: all run lclean save compile

PROG = ns
PETSCFEM_DIR   = ../..
#MYOBJS = ns.o fracstep.o nsi_tet.o
MYOBJS = ns.o nsi_tet.o nsitetles.o nsitetlesfm2.o \
		bccnsfm2.o wall.o nsvol.o walldata.o

TAGFLAGS = --include="../../src/TAGS"

# for makedepend
SRCS := $(MYOBJS:.o=.cpp)

include $(PETSCFEM_DIR)/Makefile.base

$(PROG) compile: $(PROG).bin

%.cppi: %.cpp
	g++ $(CPPFLAGS) -I /usr/local/mpich-1.2.0/include/ \
		-I/home/mstorti/PETSC/petsc-2.0.24/bmake/linux/ \
		-I/home/mstorti/PETSC/petsc-2.0.24/include/ -E $< > $@

# documentation of options with the script `odoc.pl'
odocns: ../../doc/odocns.tex ../../doc/odocnse.tex \
	../../doc/odocnsb.tex  ../../doc/odocnsw.tex 

../../doc/odocns.tex: ns.cpp force
	-chmod +w $@
	$(ODOC) -o $@ $<
	-chmod -w $@

../../doc/odocnse.tex: nsitetlesfm2.cpp force
	-chmod +w $@
	$(ODOC) -o $@ $<
	-chmod -w $@

../../doc/odocnsb.tex: bccnsfm2.cpp force
	-chmod +w $@
	$(ODOC) -o $@ $<
	-chmod -w $@

../../doc/odocnsw.tex: wall.cpp force
	-chmod +w $@
	$(ODOC) -o $@ $<
	-chmod -w $@

# $(PROG).bin: $(MYOBJS) $(LIBPETSCFEM)
# 	if [ -e $(PROG).bin ] ; then rm $(PROG).bin ; fi
# 	${CXX_CLINKER} -o $(PROG).bin $(MYOBJS) \
# 		-L../../src -lpetscfem \
#                 -L$(METIS) -L/usr/local/lib \
#                 -lmetis  -L$(NEWMAT)  -lnewmat \
# 		-lglib -L$(LIBRETTO_LIB) -libretto \
# 		-lc ${PETSC_SLES_LIB} \
# 		$(MESCHACH)/libmes.a

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

lclean: general_clean
#	cat /dev/null >makefile.d
#	-rm *.o $(OBJS) $(PROG).bin octave-core *.cppi outvector.sal makefile.d
#	$(MAKE) clean

# File containing dependencies computed by makedepend
include makefile.d
