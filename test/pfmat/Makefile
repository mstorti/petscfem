#__INSERT_LICENSE__
# $Id: Makefile,v 1.1.2.7 2002/01/03 21:47:44 mstorti Exp $ 

#s Local Targets
#w all tests
tests: pfmat pfmat2

SRCS = pfmat.cpp pfmat2.cpp 

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

CLEAN_DIRS = 
DEPEND_DIRS = $(CLEAN_DIRS)

SRC_DIR = $(PETSCFEM_DIR)/src
LIB := $(SRC_DIR)/libpetscfem.a
LIB_OBJS := graph util2 iisdmat petscmat sparse spsolve distmat iisdcr \
		iisdgraph spline seval debug texthash utils util2 fm2eperl \
		fastmat2 fstack fem getprop mainutl dofmap idmap tempfun
LIB_OBJS_LIB := $(patsubst %,$(LIB)(%.o),$(LIB_OBJS))

$(LIB)(%.o):
	$(MAKE) -C $(SRC_DIR) 'libpetscfem.a($*.o)'

#w Build executable
pfmat.bin: $(LIB_OBJS_LIB) pfmat.o
	-rm -f $@
	${CXX_CLINKER} -o $@ pfmat.o  $(LDFLAGS)

NP := $(shell procsel $(PETSCFEM_PROCTABLE))

pfmat = out=output.$(3).tmp ; rm -f $$out ;	\
	echo args: $(2) > $$out ;		\
	$(MPIRUN) -np $(1) -machinefile		\
        machi.dat ./pfmat.bin $(2) >> $$out

#w Several tests PFMat matrices
pfmat:  #pfmat.bin
#	Try all combinations of: np=1,2, iisd_subpart=1,10, mtyp=0,1
	$(call pfmat,1,ne 100 deb 0 nsl 2 sbdp 10 nmat 2 rrhs 1 mtyp 0 q 1,case1)
	$(call pfmat,2,ne 100 deb 0 nsl 2 sbdp 10 nmat 2 rrhs 1 mtyp 0 q 1,case2)
	$(call pfmat,1,ne 100 deb 0 nsl 2 sbdp 1  nmat 2 rrhs 1 mtyp 0 q 1,case3)
	$(call pfmat,2,ne 100 deb 0 nsl 2 sbdp 1  nmat 2 rrhs 1 mtyp 0 q 1,case4)
#
	$(call pfmat,1,ne 100 deb 0 nsl 2 sbdp 10 nmat 2 rrhs 1 mtyp 1 q 1,case5)
	$(call pfmat,2,ne 100 deb 0 nsl 2 sbdp 10 nmat 2 rrhs 1 mtyp 1 q 1,case6)
	$(call pfmat,1,ne 100 deb 0 nsl 2 sbdp 1  nmat 2 rrhs 1 mtyp 1 q 1,case7)
	$(call pfmat,2,ne 100 deb 0 nsl 2 sbdp 1  nmat 2 rrhs 1 mtyp 1 q 1,case8)

#w Build executable
pfmat2.bin: $(LIB_OBJS_LIB) pfmat2.o
	-rm -f $@
	${CXX_CLINKER} -o $@ pfmat2.o  $(LDFLAGS)



#s
