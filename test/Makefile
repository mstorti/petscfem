#__INSERT_LICENSE__
# $Id: Makefile,v 1.100 2003/12/08 23:24:45 mstorti Exp $ 

.PHONY: all run  compile lclean ns laplace adv nso laplace_tests force		\
        ctest runtests clean fastmat2_obj plano idmap fastvec			\
	chunk_size_process fastmat2 fastmat2a fastmat2b fastmat2c sqcav5	\
	constraint_bug fstack sector triangle advdif texthash log_vars		\
	newff burgers bug100 nonlr genload lupart distmap sparse graph sqcav	\
	pfmat  aquifer nutatcyl mmove advec nsgenload plate gasflow oscplate \
	idmap sw_abso

SRCS = $(wildcard *.cpp)

DIRS = circ fstack sector les turbchan turbchanw advdif texthash		\
		newff burgers les/stabi genload bug100 nonlr wallke		\
		wkepipe wallkechan lupart distmap sparse graph sqcav		\
		pfmat cubcav aquifer bubbly nutatcyl aquifer fastmat2		\
		nsgenload plate gasflow idmap mmove advec plano oscplate 	\
		spillway sw_abso

CLEAN_DIRS = $(DIRS)
DEPEND_DIRS = texthash distmap sparse aquifer plate fstack idmap

SWDIRS:= 
local_sw::
	chmod 755 `find . -name 'gmv*.pl'`

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
	-runtests.pl -o | tee -a make-tests-$(BOPT).log
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
	$(MAKE) BOPT=g_c++ torture_1
	$(MAKE) lclean
	$(MAKE) BOPT=O_c++ torture_1
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
runtests: compile sector oscplate fastmat2 plano idmap fastvec			\
	chunk_size_process constraint_bug fstack triangle advdif texthash	\
	log_vars newff burgers bug100 genload lupart distmap sparse graph	\
	sqcav pfmat aquifer nutatcyl mmove advec nsgenload oscplate idmap	\
	sw_abso

rp: runtests_p
runtests_p: sw_abso

#s Specific tests
#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
#w Check bug100 (don't allocate chunk_size see notes.txt)
bug100: $(ADVDIF)
	$(MAKE) -C bug100 run

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
#w Checks the constraint bug feature
constraint_bug: laplace
	$(MAKE) -C circ constraint_bug

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
#w FileStack class check
fstack:
	$(MAKE) -C fstack tests

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
	$(MAKE) -C sparse tests

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
#w Oscillating plate tests
idmap:
	$(MAKE) -C $@ tests

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>--
#w Make advective plano tests
plano:
	$(MAKE) -C $@ tests

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
#w Laplace on a sector of a circular sector
#w with periodic b.c.'s
sector: 
	$(MAKE) -C $@ tests

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>--
#w Gneric load elemeset for NS
nsgenload:
	$(MAKE) -C $@ tests

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>--
# sw absorbent bc
sw_abso:
	$(MAKE) -C sw_abso tests

#s
