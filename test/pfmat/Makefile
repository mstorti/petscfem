#__INSERT_LICENSE__
# $Id: Makefile,v 1.1.2.8 2002/01/05 23:50:26 mstorti Exp $ 
#s Local Targets
#w all tests
tests: pfmat

SRCS = pfmat.cpp pfmat2.cpp 

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

CLEAN_DIRS = 
DEPEND_DIRS = $(CLEAN_DIRS)

LIB_OBJS := graph util2 iisdmat petscmat sparse spsolve distmat iisdcr \
		iisdgraph spline seval debug texthash utils util2 fm2eperl \
		fastmat2 fstack fem getprop mainutl dofmap idmap tempfun

#w Build executable
pfmat.bin: $(LIB_OBJS_LIB) pfmat.o
	-rm -f $@
	${CXX_CLINKER} -o $@ pfmat.o  $(LDFLAGS)

NP := $(shell procsel $(PETSCFEM_PROCTABLE))

args = ne 40 deb 0 nsl 2 $(1) nmat 2 rnd 1 $(2) q 1 $(3)

pfmat = out=output.$(3).tmp ; rm -f $$out ;	\
	echo args: $(2) > $$out ;		\
	$(MPIRUN) -np $(1) -machinefile		\
        machi.dat ./pfmat.bin $(2)  >> $$out

#w Several tests PFMat matrices
pfmat:  pfmat.bin
#	Try all combinations of: np=1,2, iisd_subpart=1,10, mtyp=0,1
	$(call pfmat,1,$(call args,sbdp 4, mtyp 0, ls 0),case1)
	$(call pfmat,2,$(call args,sbdp 4, mtyp 0, ls 0),case2)
	$(call pfmat,1,$(call args,sbdp 1, mtyp 0, ls 0),case3)
	$(call pfmat,2,$(call args,sbdp 1, mtyp 0, ls 0),case4)
#		       	       	      
	$(call pfmat,1,$(call args,sbdp 4, mtyp 1, ls 0),case5)
	$(call pfmat,2,$(call args,sbdp 4, mtyp 1, ls 0),case6)
	$(call pfmat,1,$(call args,sbdp 1, mtyp 1, ls 0),case7)
	$(call pfmat,2,$(call args,sbdp 1, mtyp 1, ls 0),case8)
#
	$(call pfmat,1,$(call args,sbdp 4, mtyp 0, ls 1),case1s)
	$(call pfmat,2,$(call args,sbdp 4, mtyp 0, ls 1),case2s)
	$(call pfmat,1,$(call args,sbdp 1, mtyp 0, ls 1),case3s)
	$(call pfmat,2,$(call args,sbdp 1, mtyp 0, ls 1),case4s)
#		       	       	      
	$(call pfmat,1,$(call args,sbdp 4, mtyp 1, ls 1),case5s)
	$(call pfmat,2,$(call args,sbdp 4, mtyp 1, ls 1),case6s)
	$(call pfmat,1,$(call args,sbdp 1, mtyp 1, ls 1),case7s)
	$(call pfmat,2,$(call args,sbdp 1, mtyp 1, ls 1),case8s)

pfmat2: pfmat2.bin
	pfmat2.bin

#w Build executable
pfmat2.bin: $(LIB_OBJS_LIB) pfmat2.o
	-rm -f $@
	${CXX_CLINKER} -o $@ pfmat2.o  $(LDFLAGS)

#s
