#__INSERT_LICENSE__
# $Id: Makefile,v 1.85 2003/02/17 12:39:28 mstorti Exp $ 

.PHONY: all run  compile lclean ns laplace adv nso laplace_tests force		\
        ctest runtests clean fastmat2_obj plano idmap fastvec			\
	chunk_size_process fastmat2 fastmat2a fastmat2b fastmat2c sqcav5	\
	constraint_bug fstack sector triangle advdif texthash log_vars		\
	newff burgers bug100 nonlr genload lupart distmap sparse graph sqcav	\
	pfmat  aquifer nutatcyl mmove advec nsgenload plate gasflow oscplate

SRCS = $(wildcard *.cpp)

DIRS = circ fstack sector les turbchan turbchanw advdif texthash		\
		newff burgers les/stabi genload bug100 nonlr wallke		\
		wkepipe wallkechan lupart distmap sparse graph sqcav		\
		pfmat cubcav aquifer bubbly nutatcyl aquifer fastmat2 mmove advec \
		nsgenload plate gasflow

CLEAN_DIRS = $(DIRS)
DEPEND_DIRS = texthash distmap sparse aquifer plate fstack

run: sqcav5

include Makefile.base.test

#s Main targets
#w Makes all tests and report succes
all: tests

PETSCFEM_DIR = ..
#PFEMUP = ../..
include $(PETSCFEM_DIR)/Makefile.base

#w Builds all applications
compile:
	-rm -f *.sal
	@echo "Running tests..."
	$(MAKE) -C ../ ns laplace adv advdif

#w clean local files
local_clean::
	-rm -f ns advdif adv laplace

.PHONY: torture_1 torture_g torture_O torture

#---<*>---//---<*>---//---<*>---//---<*>---//---<*>---// 
#w Make all tests for a particular BOPT and keep log
#w and make a tar with the results
torture_1:
	$(MAKE) -k runtests 2>&1 | tee make-tests-$(BOPT).log
	-runtests.pl | tee make-tests-$(BOPT).log
	tar czvf ../make-tests-$(BOPT).tar.gz .

#w torture for BOPT=g_c++
torture_g:
	$(MAKE) BOPT=g_c++ torture_1

#w torture for BOPT=O_c++
torture_O:
	$(MAKE) BOPT=O_c++ torture_1

#w torture for BOPT=g_c++ and BOPT=g_c++,
#w store results in tars and logs.
torture:
	$(MAKE) lclean
	$(MAKE) BOPT=O_c++ torture_1
	$(MAKE) lclean
	$(MAKE) BOPT=g_c++ torture_1
	$(MAKE) lclean

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
# All tests

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
#w Verifies results
tests: runtests
	@echo "Verifying results..."
	@runtests.pl
	@echo "Done."

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
#w Runs the tests
runtests: compile lap_per sector oscplate fastmat2 plano idmap fastvec		\
	chunk_size_process constraint_bug fstack triangle advdif texthash	\
	log_vars newff burgers bug100 genload lupart distmap sparse graph	\
	sqcav pfmat aquifer nutatcyl mmove advec nsgenload oscplate

runtests_p:   

#s Specific tests
#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
#w Lapace on a square with periodic b.c.'s
lap_per: lap_per.dat laplace
	-$(MPIRUN) -np 1 $(LAPLACE) -case lap_per.dat >/dev/null

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
#w Check bug100 (don't allocate chunk_size see notes.txt)
bug100: $(ADVDIF)
	$(MAKE) -C bug100 run

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
#w Laplace on a sector of a circular sector
#w with periodic b.c.'s
sector: 
	$(MAKE) -C sector sector

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
.PHONY: $(plano_list)
plano_list = plano_cons plano_cons_auto plano_local plano_auto plano_local_nm \
		plano_local_weak plano_fm2_weak 

#w Test for shallow water
plano: adv $(plano_list)

plano.depl: force

ADVO = $(PETSCFEM_DIR)/../petscfem-beta-1.14/applications/advective/adv.bin

export EPERLFLAGS

ifneq ($(RH60),) 

plano_cons:       EPERLFLAGS = "-d consistent_supg=1"
plano_cons_auto:  EPERLFLAGS = "-d consistent_supg_auto=1"
plano_local:      EPERLFLAGS = "-d local_time_step=1"
plano_auto:       EPERLFLAGS = "-d auto_time_step=1"
plano_local_nm:   EPERLFLAGS = "-d local_time_step_nm=1"
plano_local_weak: EPERLFLAGS = "-d local_time_step_weak=1"
plano_fm2_weak:   EPERLFLAGS = "-d fm2_weak=1"

plano_%: 
	$(MAKE) EPERLFLAGS=$(EPERLFLAGS) plano.depl
	-$(MPIRUN) -np 1 $(ADV) -case plano.depl &> $@.sal
else

plano_cons:
	$(MAKE) EPERLFLAGS="-d consistent_supg=1" plano.depl
	-$(MPIRUN) -np 1 $(ADV) -case plano.depl &> $@.sal

plano_cons_auto:
	$(MAKE) EPERLFLAGS="-d consistent_supg_auto=1" plano.depl
	-$(MPIRUN) -np 1 $(ADV) -case plano.depl &> $@.sal

plano_local:      
	$(MAKE) EPERLFLAGS="-d local_time_step=1" plano.depl
	-$(MPIRUN) -np 1 $(ADV) -case plano.depl &> $@.sal

plano_auto:       
	$(MAKE) EPERLFLAGS="-d auto_time_step=1" plano.depl
	-$(MPIRUN) -np 1 $(ADV) -case plano.depl &> $@.sal

plano_local_nm:
	$(MAKE) EPERLFLAGS="-d local_time_step_nm=1" plano.depl
	-$(MPIRUN) -np 1 $(ADV) -case plano.depl &> $@.sal

plano_local_weak: 
	$(MAKE) EPERLFLAGS="-d local_time_step_weak=1" plano.depl
	-$(MPIRUN) -np 1 $(ADV) -case plano.depl &> $@.sal

plano_fm2_weak:   
	$(MAKE) EPERLFLAGS="-d fm2_weak=1" plano.depl
	-$(MPIRUN) -np 1 $(ADV) -case plano.depl &> $@.sal

endif


#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
#w Idmap class
idmap: tidmap.bin
	tidmap.bin >tidmap.sal

../src/%.o: force
	$(MAKE) -C ../src $@

# tidmap.bin: tidmap.o ../src/idmap.o ../src/utils.o 
# #	g++ -o tidmap.bin   -L$(NEWMAT) -lnewmat
# 	${CXX_CLINKER} -o $@ tidmap.o ../src/idmap.o ../src/utils.o \
# 	../src/fastmat2.o ../src/fm2eperl.o \
# 		-lc ${PETSC_SLES_LIB} -L$(NEWMAT) -lnewmat

tidmap.bin: tidmap.o $(LIBPETSCFEM)
#	g++ -o tidmap.bin   -L$(NEWMAT) -lnewmat
	${CXX_CLINKER} -o $@ tidmap.o $(LIBPETSCFEM) \
		-lc ${PETSC_SLES_LIB} -L$(NEWMAT) -lnewmat  $(MESCHACH)/libmes.a



#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
#w Checks the constraint bug feature
constraint_bug: laplace
	$(MAKE) -C circ constraint_bug

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
#w 5x5 Square cavity for NS
NSO = $(PETSCFEM_DIR)/../petscfem-2.16g/applications/ns/ns.bin
sqcav5_dat = sqcav5.nod.tmp sqcav5.con.tmp sqcav5.fixa.tmp
sqcav5: ns sqcav5_dat
	-$(MPIRUN) -np 3 $(NSO) -case sqcav5.dat | tee sqcav5.out

sqcav5_dat: mksqcav5.m
	$(OCTAVE)  mksqcav5.m >/dev/null

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
#w Tests processing of chunks
chunk_size_process:
	$(MAKE) -C sector chunk_size_process

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
#w Tests laplace with triangles
triangle:
	$(MAKE) -C sector triangle

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
#w FileStack class check
fstack:
	$(MAKE) -C fstack fstack

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
#w Texthash class check
texthash:
	$(MAKE) -C texthash tests

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
# Tests the new iterator class
iterator: iterator.cpp
	g++ -o $@ $<

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
#w Tests for the advdif module
advdif:
	$(MAKE) -C advdif tests
	$(MAKE) -C turbchan turbchan

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
#w Logarithmic variables
log_vars:
	$(MAKE) -C turbchanw tests

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
#s FastMat2 tests
#w All FastMat2 tests
fastmat2: 
	$(MAKE) -C fastmat2 tests

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>--
#w New flux functions
newff:
	$(MAKE) -C newff tests

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>--
#w Generic load element for advdif
genload:
	$(MAKE) -C genload tests

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>--
#w IISD solver and alternative partitionings
lupart:
	$(MAKE) -C lupart tests

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>--
#w Specific tests for the DistMap class/template
distmap:
	$(MAKE) -C distmap tests

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>--
#w Burgers with new flux functions
burgers:
	$(MAKE) -C burgers tests

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>--
#w Sparse matrices
sparse:
	$(MAKE) -C sparse sparse

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>--
#w Graph class
graph:
	$(MAKE) -C graph tests

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>--
#w Square cavity tests
sqcav:
	$(MAKE) -C sqcav tests

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>--
#w Tests on the PFMat class
pfmat:
	$(MAKE) -C pfmat tests

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>--
#w Tests for the Hydrology related elemsets
aquifer:
	$(MAKE) -C aquifer tests

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>--
#w Navier Stokes test for the nutating cylinder
nutatcyl:
	$(MAKE) -C nutatcyl tests

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>--
#w Mesh momement tests
mmove:
	$(MAKE) -C mmove tests

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>--
#w Linear advection 
advec:
	$(MAKE) -C $@ tests

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>--
#w Oscillating plate tests
oscplate:
	$(MAKE) -C $@ tests

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>--
#w Gneric load elemeset for NS
nsgenload:
	$(MAKE) -C $@ tests

#s
