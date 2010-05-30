#__INSERT_LICENSE__

.PHONY: all run  compile lclean ns laplace adv nso laplace_tests force		\
        ctest runtests clean fastmat2_obj plano idmap				\
	chunk_size_process fastmat2 fastmat2a fastmat2b fastmat2c sqcav5	\
	constraint_bug fstack sector triangle advdif texthash log_vars		\
	newff burgers bug100 nonlr genload lupart distmap sparse graph sqcav	\
	pfmat  aquifer nutatcyl mmove advec nsgenload plate gasflow-conv 	\
	oscplate idmap sw_abso cubcav condwall project embgath alebcconv 	\
	ductadv

SRCS = $(wildcard *.cpp)

DIRS = circ fstack sector les turbchan turbchanw advdif texthash	\
		newff burgers les/stabi genload bug100 nonlr wallke	\
		wkepipe wallkechan lupart distmap sparse graph sqcav	\
		pfmat cubcav aquifer bubbly nutatcyl aquifer fastmat2	\
		nsgenload plate gasflow-conv idmap mmove advec plano	\
		oscplate spillway sw_abso gfplate condwall project	\
		embgath flow-reversal alebcconv ductadv

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
	mpdboot
	$(MAKE) -k runtests 2>&1 | tee make-tests-$(BOPT).log
#
	echo "-------------------------" | tee -a  make-tests-$(BOPT).log
	echo "FULL TEST REPORT: " | tee -a make-tests-$(BOPT).log
	-./runtests.pl -o | tee -a make-tests-$(BOPT).log
#
	echo "-------------------------" | tee -a  make-tests-$(BOPT).log
	echo "COULDN'T OPEN: " | tee -a make-tests-$(BOPT).log
	-./runtests.pl -o | grep "^can't" | tee -a make-tests-$(BOPT).log
#
	echo "-------------------------" | tee -a  make-tests-$(BOPT).log
	echo "NOT OK TESTS: " | tee -a make-tests-$(BOPT).log
	-./runtests.pl -n | tee -a make-tests-$(BOPT).log
#
	tar czvf ../make-tests-$(BOPT).tar.gz -C .. test
	echo "-------------------------" | tee -a make-tests-$(BOPT).log
	echo "ABRIDGED: " | tee -a  make-tests-$(BOPT).log
	$(MAKE) abridged &> /dev/null
	cat abridged.tmp >> make-tests-$(BOPT).log

#w report:
report:
	echo "-------------------------" | tee -a  make-tests-$(BOPT).log
	echo "FULL TEST REPORT: " | tee -a make-tests-$(BOPT).log
	-./runtests.pl -o | tee -a make-tests-$(BOPT).log
#
	echo "-------------------------" | tee -a  make-tests-$(BOPT).log
	echo "COULDN'T OPEN: " | tee -a make-tests-$(BOPT).log
	-./runtests.pl -o | grep "^can't" | tee -a make-tests-$(BOPT).log
#
	echo "-------------------------" | tee -a  make-tests-$(BOPT).log
	echo "NOT OK TESTS: " | tee -a make-tests-$(BOPT).log
	-./runtests.pl -n | tee -a make-tests-$(BOPT).log

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
	mkstepsdir.pl -p petscfem-runtests
	mv -f make-tests* petscfem-runtests
	mv -f ../make-tests-*.tar.gz petscfem-runtests

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
# All tests

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
#w Verifies results
tests: runtests
	@echo "Verifying results..."
	@./runtests.pl
	@echo "Done."

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
#w Runs the tests
runtests: compile sector oscplate fastmat2 plano idmap constraint_bug	\
	fstack advdif texthash log_vars newff burgers bug100 genload	\
	lupart distmap sparse graph sqcav pfmat aquifer nutatcyl mmove	\
	advec nsgenload cubcav project embgath alebcconv ductadv 	\
	gasflow-conv

rp: runtests_p
runtests_p: pfmat aquifer nutatcyl mmove	\
	advec nsgenload cubcav project embgath alebcconv ductadv 	\
	gasflow-conv

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
	$(MAKE) -C flow-reversal tests

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>--
# sw absorbent bc
sw_abso:
	$(MAKE) -C sw_abso tests

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>--
# Tests in cubic (3D) cavity
cubcav:
	$(MAKE) -C cubcav tests

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>--
# Projections of data between teo meshes
project:
	$(MAKE) -C ../tools/project project4.bin
	$(MAKE) -C project tests

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>--
# Embedded gatherer test
embgath:
	$(MAKE) -C embgath tests

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>--
# Invariance of ALE with bcconv (ALE corrected version)
alebcconv:
	$(MAKE) -C alebcconv tests

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>--
# Check strong temporal conservation with advdif
ductadv:
	$(MAKE) -C ductadv tests

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>--
# Gasflow tests
gasflow-conv:
	$(MAKE) -C $@ tests

# makes an abridged acount of the tests
abridged2:
	@ ( echo "NOT OK-------" ;				\
	./runtests.pl -n | grep "^Testing" |			\
		perl -pe 's/^.*on file "(.*)".*$$/$$1/' ;	\
	echo "CAN'T OPEN-------" ;				\
	./runtests.pl -o | grep "^can't open" |			\
		perl -pe "s/^can't open //" ; )			\
		| tee abridged.tmp

abridged:
	cat /dev/null > abridged.tmp
	./runtests.pl -n | grep "^Testing" |			\
		perl -pe 's/^.*on file "(.*)".*$$/$$1/' > tempo.txt
	echo -n "NOT OK------- -> " >> abridged.tmp
	wc -l tempo.txt | cut -f 1 -d ' '>> abridged.tmp
	cat tempo.txt >> abridged.tmp
#
	./runtests.pl -o | grep "^can't open" |			\
		perl -pe "s/^can't open //" > tempo.txt
	echo -n "CAN'T OPEN------- -> " >> abridged.tmp
	wc -l tempo.txt | cut -f 1 -d ' '>> abridged.tmp
	cat tempo.txt >> abridged.tmp
#
	rm -f tempo.txt

#s
