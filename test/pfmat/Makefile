#__INSERT_LICENSE__
#s Local Targets

default: pfmat7

#w all tests
# tests: pfmat pfmat2
tests: pfmat

.PHONY: pfmat pfmat2 force iisdmat_petsc iisdmat_slu petscmat \
			spd_petsc spd_slu

SRCS = pfmat.cpp pfmat2.cpp 

LIB_OBJS := graph util2 iisdmat petscmat sparse spsolve distmat iisdcr \
		iisdgraph spline seval debug texthash utils util2 fm2eperl \
		fastmat2 fstack fem getprop mainutl dofmap idmap tempfun \
		spdirect

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

LIB := $(SRC_DIR)/libpetscfem$(osfx).a
LIB_OBJS_LIB := $(patsubst %,$(LIB)(%.o),$(LIB_OBJS))

$(LIB)(%.o):
	$(MAKE) -C $(SRC_DIR) 'libpetscfem$(osfx).a($*.o)'

CLEAN_DIRS = 
DEPEND_DIRS = $(CLEAN_DIRS)

#w Build executable
pfmat.bin: $(LIB_OBJS_LIB) pfmat.o
	-rm -f $@
	${CXX_CLINKER} -o $@ pfmat.o  $(LDFLAGS_NO_PROG_LIB)

NP := $(shell procsel $(PETSCFEM_PROCTABLE))

args = ne 20 deb 0 nsl 2 $(1) nmat 2 rnd 1 $(2) q 1 $(3)

# MCH := -machinefile machi.dat

pfmat = out=output.$(3).tmp ; rm -f $$out ;	\
	echo   args: $(2) > $$out ;		\
	$(MPIEXEC) $(MCH) -n $(1)	\
	./pfmat.bin $(2) </dev/null >> $$out

pfmat: pfmat.bin iisdmat_petsc petscmat spd_petsc spd_slu

#w IISDMat/PETSc  matrix tests
iisdmat_petsc:
	$(call pfmat,1,$(call args,sbdp 4, mtyp 0),case1iip)
	$(call pfmat,2,$(call args,sbdp 4, mtyp 0),case2iip)
	$(call pfmat,1,$(call args,sbdp 1, mtyp 0),case3iip)
	$(call pfmat,2,$(call args,sbdp 1, mtyp 0),case4iip)

#w IISDMat/SuperLU  matrix tests
iisdmat_slu:
	$(call pfmat,1,$(call args,sbdp 4, mtyp 1),case1iis)
	$(call pfmat,2,$(call args,sbdp 4, mtyp 1),case2iis)
	$(call pfmat,1,$(call args,sbdp 1, mtyp 1),case3iis)
	$(call pfmat,2,$(call args,sbdp 1, mtyp 1),case4iis)

#w PETScMat  matrix tests
petscmat:
	$(call pfmat,1,$(call args,sbdp 1, mtyp 2),case1p)
	$(call pfmat,2,$(call args,sbdp 1, mtyp 2),case2p)

#w SparseDirect/PETSc  matrix tests
spd_petsc:
	$(call pfmat,1,$(call args,sbdp 1, mtyp 3),case1sdp)

#w SparseDirect/Superlu matrix tests
spd_slu:
	$(call pfmat,1,$(call args,sbdp 1, mtyp 4),case1sds)

pp:
#	Try all combinations of: np=1,2, iisd_subpart=1,4, mtyp=0,1
	for np in 1 2 ; do					\
		for sbdp in 1 4 ; do				\
			for mtyp in 0 1 2 3 4 ; do		\
		if [ $$np -eq 1 -o $$mtyp -le 2  ] ; then	\
	out=output.np$$np.sbpd$$sbdp.mtyp$$mtyp.tmp ;		\
        rm -f $$out ;						\
	args="ne 20 deb 0 nsl 2 sbdp $$sbdp  nmat 2 rnd		\
		1 mtyp $$mtyp q 1" ;				\
	echo  args: $$args > $$out ;				\
	$(MPIEXEC) -n $$np $(MCH)				\
		./pfmat.bin $$args </dev/null  >> $$out  ;	\
		fi ;						\
			done ;					\
		done ;						\
	done 

args2 = ne 100 deb 0 nsl 2 $(1) nmat 2 nprof 2 $(2) f 0.9

pfmat2 = out=output2.$(3).tmp ; rm -f $$out ;		\
	echo args: $(2) > $$out ;			\
	$(MPIEXEC) $(MCH) -n $(1)  ./pfmat2.bin $(2)	\
		</dev/null  >> $$out

pfmat2: pfmat2.bin force
	$(call pfmat2,1,$(call args2,sbdp 4, mtyp 0, ls 0),case1)
	$(call pfmat2,2,$(call args2,sbdp 4, mtyp 0, ls 0),case2)
	$(call pfmat2,1,$(call args2,sbdp 1, mtyp 0, ls 0),case3)
	$(call pfmat2,2,$(call args2,sbdp 1, mtyp 0, ls 0),case4)
#
	$(call pfmat2,1,$(call args2,sbdp 4, mtyp 0, ls 1),case1s)
	$(call pfmat2,2,$(call args2,sbdp 4, mtyp 0, ls 1),case2s)
	$(call pfmat2,1,$(call args2,sbdp 1, mtyp 0, ls 1),case3s)
	$(call pfmat2,2,$(call args2,sbdp 1, mtyp 0, ls 1),case4s)

#w Build executable
pfmat2.bin: $(LIB_OBJS_LIB) pfmat2.o
	-rm -f $@
	${CXX_CLINKER} -o $@ pfmat2.o  $(LDFLAGS_NO_PROG_LIB)

#w Build executable
pfmat4.bin: $(LIB_OBJS_LIB) pfmat4.o
	-rm -f $@
	${CXX_CLINKER} -o $@ pfmat4.o  $(LDFLAGS_NO_PROG_LIB)

pfmat5.bin: pfmat5.o
	$(MAKE) -C ../../src compile
	-rm -f $@
	${CXX_CLINKER} -o $@ pfmat5.o  $(LDFLAGS_NO_PROG_LIB)
	-rm -f pfmat5.o

pfmat7.bin: pfmat7.o
	$(MAKE) -C ../../src compile
	-rm -f $@
	${CXX_CLINKER} -o $@ pfmat7.o  $(LDFLAGS_NO_PROG_LIB)
	-rm -f pfmat5.o

sw_nn_preco.bin: sw_nn_preco.o
	$(MAKE) -C ../../src compile
	-rm -f $@
	${CXX_CLINKER} -o $@ sw_nn_preco.o  $(LDFLAGS_NO_PROG_LIB)
	-rm -f sw_nn_preco.o

try.bin: try.o
	-rm -f $@
	${CXX_CLINKER} -o $@ try.o  $(LDFLAGS_NO_PROG_LIB)


debug.o pfmat2.o: CPPFLAGS += -O0

# CPPFLAGS += -I/usr/include/g++-3

pfmat3.bin: $(LIB_OBJS_LIB) pfmat3.o hashf.o
	-rm -f $@
	${CXX_CLINKER} -o $@ pfmat3.o hashf.o $(LDFLAGS_NO_PROG_LIB)

pfmat7: pfmat7.bin
	$(MPIEXEC) -n 1 ./pfmat7.bin $(call in_backg,$@)

#s
