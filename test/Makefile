#__INSERT_LICENSE__
# $Id: Makefile,v 1.25 2001/05/27 17:15:23 mstorti Exp $ 

.PHONY: all run  compile lclean ns laplace adv nso laplace_tests force \
        ctest runtests clean fastmat2_obj plano idmap fastvec \
	chunk_size_process fastmat2 fastmat2a fastmat2b fastmat2c sqcav5 \
	constraint_bug fstack sector triangle advdif texthash log_vars \
	newff burgers bug100

run: sqcav5

include Makefile.base.test

#s Main targets
#w Makes all tests and report succes
all: tests

PETSCFEM_DIR = ..
#PFEMUP = ../..
include $(PETSCFEM_DIR)/Makefile.base

SRCS = $(wildcard *.cpp)

DIRS = circ fstack sector les turbchan turbchanw advdif texthash \
		newff burgers les/stabi genload bug100
CLEAN_DIRS = $(DIRS)
DEPEND_DIRS = $(DIRS)

ctest: $(DATAFILE) force
	$(MAKE) $(NS)
	-$(MPIRUN) -np 1 $(NS) -case oscplate.depl

#w Builds all applications
compile:
	-rm *.sal
	@echo "Running tests..."
	$(MAKE) -C ../ ns laplace adv advdif

#w clean local files
local_clean::
	-rm -f ns advdif adv laplace

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
runtests: compile lap_per sector oscplate_all fastmat2 plano \
	idmap fastvec chunk_size_process constraint_bug fstack triangle \
	advdif texthash log_vars newff burgers bug100
runtests_p: 	idmap fastvec chunk_size_process constraint_bug fstack triangle \
	advdif texthash log_vars newff burgers

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
#w All oscillating plate tests
#w saco oscplate2b 
oscplate_all: oscplate3 oscplate2 oscplate2b oscplate4 \
	oscplate4cos oscplate4_30deg oscplate

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
#w Oscillating plate with piecewise constant function
oscplate2: ns
	-rm tmp_file*.tmp oscsome2.sal
	eperl -P -d piecewise_fun=1  oscplate.epl >oscplate2.depl
	-$(MPIRUN) -np 1 $(NS) -case oscplate2.depl >/dev/null
	octave oscplate2.m &> oscplate2o.sal

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
#w Oscillating plate with sinusoidal oscillation
oscplate: ns
	-rm tmp_file*.tmp oscplate1.sal
	eperl -P -d Nperiods=20 oscplate.epl >oscplate.depl
	-$(MPIRUN) -np 1 $(NS) -case oscplate.depl >/dev/null
	mv oscplate.sal oscplate1.sal

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
#w Oscillating plate with sinusoidal oscillation
#w Verifies quadratic convergence of the
#w Crank Nicholson method. Runs a case with N=16,
#w 32 and 128 time steps. The solution for N=128 is
#w considered the exact, and it is checked in Octave that
#w norm(u32-u128)/norm(u16-u128) < 0.30 (it should be <0.25). 
oscplate3: ns
	-rm tmp_file*.tmp oscsome3.sal oscsome.sal
	-chmod +w oscplate3.depl
#
	eperl -P -d alpha=0.5 -d N=16 oscplate.epl >oscplate3.depl
	-$(MPIRUN) -np 1 $(NS) -case oscplate3.depl  >/dev/null
	mv oscplate.sal oscplate_16.sal
#
	eperl -P -d alpha=0.5 -d N=32 oscplate.epl >oscplate3.depl
	-$(MPIRUN) -np 1 $(NS) -case oscplate3.depl  >/dev/null
	mv oscplate.sal oscplate_32.sal
#
	eperl -P -d alpha=0.5 -d N=128 oscplate.epl >oscplate3.depl
	-$(MPIRUN) -np 1 $(NS) -case oscplate3.depl  >/dev/null
	mv oscplate.sal oscplate_128.sal
#
	octave oscplate3.m &> oscplate3o.sal

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
#w Oscillating plate with piecewise constant function
oscplate2b: ns
	-rm tmp_file*.tmp oscsome2b.sal oscsome.sal
	eperl -P -d spline_fun=1  oscplate.epl >oscplate2b.depl
	-$(MPIRUN) -np 1 $(NS) -case oscplate2b.depl >/dev/null

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
#w Oscillating plate with piecewise constant periodic function
#w `sin(t/T)' case
oscplate4: ns
	-rm tmp_file*.tmp oscsome4.sal oscsome.sal
	eperl -P -d spline_periodic_fun=1 oscplate.epl >oscplate4.depl
	-$(MPIRUN) -np 1 $(NS) -case oscplate4.depl >/dev/null
	mv oscsome.sal oscsome4.sal
	octave oscplate4.m &> oscplate4o.sal

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
#w `cos(2*pi*t/T)' case
oscplate4cos: ns
	-rm tmp_file*.tmp oscsome4c.sal oscsome.sal
	eperl -P -d spline_periodic_fun=1  -d cos_case=1 \
			oscplate.epl >oscplate4c.depl
	-$(MPIRUN) -np 1 $(NS) -case oscplate4c.depl >/dev/null
	mv oscsome.sal oscsome4c.sal
	octave oscplate4c.m &> oscplate4co.sal

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
#w `sin(2*pi*t/T-pi/6)' case
oscplate4_30deg: ns
	-rm tmp_file*.tmp oscsome4_30deg.sal oscsome.sal
	eperl -P -d spline_periodic_fun=1 -d cos_case=0 -d sin30deg_case=1 \
			oscplate.epl >oscplate4_30deg.depl
	-$(MPIRUN) -np 1 $(NS) -case oscplate4_30deg.depl >/dev/null
	mv oscsome.sal oscsome4_30deg.sal
	octave oscplate4_30deg.m &> oscplate4_30dego.sal

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
#w FAST Library
fastvec: tfastvec.bin
	tfastvec.bin > tfastvec.sal

tfastvec.bin: tfastvec.o
	g++ -o tfastvec.bin tfastvec.o 

tfastvec.o: tfastvec.cpp ../src/fastlib.h

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
		-lc ${PETSC_SLES_LIB} -L$(NEWMAT) -lnewmat


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
	octave -q  mksqcav5.m

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
	$(MAKE) -C texthash texthash

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
fastmat2: fastmat2a fastmat2b fastmat2c

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>------
#w FastMat2 test a 
fastmat2a: testfm2a.bin
	testfm2a.bin > fastmat2a.sal

testfm2a.bin: testfm2a.o libpetscfem
	g++ $(PROF_FLAGS) -o testfm2a.bin testfm2a.o $(LIBPETSCFEM) \
		-L$(NEWMAT) -lnewmat ${PETSC_SLES_LIB}
# testfm2a.bin: testfm2a.o
# 	-rm $@
# 	${CXX_CLINKER} $(PROF_FLAGS) -o $@ testfm2a.o $(LDFLAGS)

daxpy.o: daxpy.f
	g77 -O3 -c daxpy.f

testfm2d: #testfm2d.bin
	./testfm2d.bin


testfm2d.bin: testfm2d.o daxpy.o libpetscfem
	g++ $(PROF_FLAGS) -o testfm2d.bin testfm2d.o daxpy.o $(LIBPETSCFEM) \
		-L$(NEWMAT) -lnewmat ${PETSC_SLES_LIB}

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>--
#w FastMat2 test b
fastmat2b: testfm2.bin
	testfm2.bin > fastmat2.sal

testfm2.bin: testfm2.o libpetscfem
	${CXX_CLINKER} $(PROF_FLAGS) -o testfm2.bin testfm2.o \
		$(LIBPETSCFEM) ${PETSC_SLES_LIB} -L$(NEWMAT) -lnewmat

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>--
#w FastMat2 test c 
fastmat2c: testfm2c.bin
	testfm2c.bin > fastmat2c.sal

testfm2c.bin: testfm2c.o libpetscfem
	g++ $(PROF_FLAGS) -o testfm2c.bin testfm2c.o $(LIBPETSCFEM) \
		 -L$(NEWMAT) -lnewmat ${PETSC_SLES_LIB}

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>--
#w New flux functions
newff:
	$(MAKE) -C newff tests

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>--
#w Burgers with new flux functions
burgers:
	$(MAKE) -C burgers tests


#s
