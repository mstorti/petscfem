# $Id: Makefile,v 1.1.1.1 2000/12/28 12:54:43 mstorti Exp $ 

.PHONY: all run  compile lclean ns laplace adv nso laplace_tests force \
        ctest runtests clean fastmat2_obj plano idmap fastvec \
	chunk_size_process fastmat2 fastmat2a fastmat2b fastmat2c sqcav5 \
	constraint_bug fstack sector triangle

include Makefile.base.test
run: sqcav5

all: tests

PETSCFEM_DIR = ..
#PFEMUP = ../..
include $(PETSCFEM_DIR)/Makefile.base

SRCS = testfm2.cpp  testfm2a.cpp testfm2d.cpp  tfastvec.cpp  tidmap.cpp \
                 tryme.cpp tryme1.cpp tryme2.cpp fstackchk.cpp

CLEAN_DIRS = circ fstack sector les turbchan advdif
DEPEND_DIRS = $(CLEAN_DIRS)

ctest: $(DATAFILE) force
	$(MAKE) $(NS)
	-$(MPIRUN) -np 1 $(NS) -case oscplate.depl

compile:
	-rm *.sal
	@echo "Running tests..."
	$(MAKE) -C ../ ns
	$(MAKE) -C ../ laplace
	$(MAKE) -C ../ adv

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
# All tests

runtests: compile lap_per sector oscplate_all fastmat2 plano \
	idmap fastvec chunk_size_process constraint_bug fstack triangle

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
LIBPFM = ../src/libpetscfem.a

$(LIBPFM): force
	$(MAKE) -C ../src compile

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
lap_per: lap_per.dat laplace
	-$(MPIRUN) -np 1 $(LAPLACE) -case lap_per.dat >/dev/null

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
sector: 
	$(MAKE) -C sector sector

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
# All oscillating plate tests
# saco oscplate2b 
oscplate_all: oscplate3 oscplate2 oscplate2b oscplate4 \
	oscplate4cos oscplate4_30deg oscplate

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
# Oscillating plate with piecewise constant function
oscplate2: ns
	-rm tmp_file*.tmp oscsome2.sal
	eperl -P -d piecewise_fun=1  oscplate.epl >oscplate2.depl
	-$(MPIRUN) -np 1 $(NS) -case oscplate2.depl >/dev/null
	octave oscplate2.m &> oscplate2o.sal

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
# Oscillating plate with sinusoidal oscillation
oscplate: ns
	-rm tmp_file*.tmp oscplate1.sal
	eperl -P -d Nperiods=20 oscplate.epl >oscplate.depl
	-$(MPIRUN) -np 1 $(NS) -case oscplate.depl >/dev/null
	mv oscplate.sal oscplate1.sal

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
# Oscillating plate with sinusoidal oscillation
# Verifies quadratic convergence of the
# Crank Nicholson method. Runs a case with N=16,
# 32 and 128 time steps. The solution for N=128 is
# considered the exact, and it is checked in Octave that
# norm(u32-u128)/norm(u16-u128) < 0.30 (it should be <0.25). 
#
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
# Oscillating plate with piecewise constant function
oscplate2b: ns
	-rm tmp_file*.tmp oscsome2b.sal oscsome.sal
	eperl -P -d spline_fun=1  oscplate.epl >oscplate2b.depl
	-$(MPIRUN) -np 1 $(NS) -case oscplate2b.depl >/dev/null

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
# Oscillating plate with piecewise constant periodic function
# `sin(t/T)' case
oscplate4: ns
	-rm tmp_file*.tmp oscsome4.sal oscsome.sal
	eperl -P -d spline_periodic_fun=1 oscplate.epl >oscplate4.depl
	-$(MPIRUN) -np 1 $(NS) -case oscplate4.depl >/dev/null
	mv oscsome.sal oscsome4.sal
	octave oscplate4.m &> oscplate4o.sal

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
# `cos(2*pi*t/T)' case
oscplate4cos: ns
	-rm tmp_file*.tmp oscsome4c.sal oscsome.sal
	eperl -P -d spline_periodic_fun=1  -d cos_case=1 \
			oscplate.epl >oscplate4c.depl
	-$(MPIRUN) -np 1 $(NS) -case oscplate4c.depl >/dev/null
	mv oscsome.sal oscsome4c.sal
	octave oscplate4c.m &> oscplate4co.sal

# #----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
# # `sin(2*pi*t/T)' case
# oscplate4cos: ns
# 	-rm tmp_file*.tmp oscsome4c.sal oscsome.sal
# 	eperl -P -d spline_periodic_fun=1  -d cos_case=1 \
# 			oscplate.epl >oscplate4c.depl
# 	-$(MPIRUN) -np 1 $(NS) -case oscplate4c.depl >/dev/null
# 	mv oscsome.sal oscsome4c.sal
# 	octave oscplate4c.m &> oscplate4co.sal

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
# `sin(2*pi*t/T-pi/6)' case
oscplate4_30deg: ns
	-rm tmp_file*.tmp oscsome4_30deg.sal oscsome.sal
	eperl -P -d spline_periodic_fun=1 -d cos_case=0 -d sin30deg_case=1 \
			oscplate.epl >oscplate4_30deg.depl
	-$(MPIRUN) -np 1 $(NS) -case oscplate4_30deg.depl >/dev/null
	mv oscsome.sal oscsome4_30deg.sal
	octave oscplate4_30deg.m &> oscplate4_30dego.sal

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
fastmat2: fastmat2a fastmat2b fastmat2c

fastmat2a: testfm2a.bin
	testfm2a.bin > fastmat2a.sal

testfm2a.bin: testfm2a.o $(LIBPFM)
	g++ $(PROF_FLAGS) -o testfm2a.bin testfm2a.o $(LIBPFM) \
		-L$(NEWMAT) -lnewmat

daxpy.o: daxpy.f
	g77 -O3 -c daxpy.f

testfm2d: #testfm2d.bin
	./testfm2d.bin


testfm2d.bin: testfm2d.o daxpy.o $(LIBPFM)
	g++ $(PROF_FLAGS) -o testfm2d.bin testfm2d.o daxpy.o $(LIBPFM) \
		-L$(NEWMAT) -lnewmat 

fastmat2b: testfm2.bin
	testfm2.bin > fastmat2.sal

testfm2.bin: testfm2.o $(LIBPFM)
	${CXX_CLINKER} $(PROF_FLAGS) -o testfm2.bin testfm2.o \
		$(LIBPFM) ${PETSC_SLES_LIB} -L$(NEWMAT) -lnewmat

fastmat2c: testfm2c.bin
	testfm2c.bin > fastmat2c.sal

testfm2c.bin: testfm2c.o $(LIBPFM)
	g++ $(PROF_FLAGS) -o testfm2c.bin testfm2c.o $(LIBPFM) \
		$(LIBPFM) -L$(NEWMAT) -lnewmat

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
# FAST Library
fastvec: tfastvec.bin
	tfastvec.bin > tfastvec.sal

tfastvec.bin: tfastvec.o
	g++ -o tfastvec.bin tfastvec.o 

tfastvec.o: tfastvec.cpp ../src/fastlib.h

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
# Test for shallow water
plano: adv
	if [ -e plano.depl ] ; then chmod +w plano.depl ; fi
	eperl -P -d consistent_supg=1  plano.epl >plano.depl
	-$(MPIRUN) -np 1 $(ADV) -case plano.depl &> plano_cons.sal
	eperl -P -d consistent_supg_auto=1  plano.epl >plano.depl
	-$(MPIRUN) -np 1 $(ADV) -case plano.depl &> plano_cons_auto.sal
	eperl -P -d local_time_step=1  plano.epl >plano.depl
	-$(MPIRUN) -np 1 $(ADV) -case plano.depl &> plano_local.sal
	eperl -P -d auto_time_step=1  plano.epl >plano.depl
	-$(MPIRUN) -np 1 $(ADV) -case plano.depl &> plano_auto.sal
	eperl -P -d local_time_step_nm=1 plano.epl >plano.depl
	-$(MPIRUN) -np 1 $(ADV) -case plano.depl &> plano_local_nm.sal
	eperl -P -d local_time_step_weak=1  plano.epl >plano.depl
	-$(MPIRUN) -np 1 $(ADV) -case plano.depl  &> plano_local_weak.sal
	eperl -P -d fm2_weak=1  plano.epl >plano.depl
	-$(MPIRUN) -np 1 $(ADV) -case plano.depl  &> plano_fm2_weak.sal

planoo:
	eperl -P -d fm2_weak=1  plano.epl >plano.depl
	-$(MPIRUN) -np 1 $(ADV) -case plano.depl  &> plano_fm2_weak.sal

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
# Idmap
idmap: tidmap.bin
	tidmap.bin >tidmap.sal

../src/%.o: force
	$(MAKE) -C ../src $@

tidmap.bin: tidmap.o ../src/idmap.o ../src/utils.o 
#	g++ -o tidmap.bin   -L$(NEWMAT) -lnewmat
	${CXX_CLINKER} -o $@ tidmap.o ../src/idmap.o ../src/utils.o \
	../src/fastmat2.o ../src/fm2eperl.o \
		-lc ${PETSC_SLES_LIB} -L$(NEWMAT) -lnewmat

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
# Checks the constraint bug feature
constraint_bug: laplace
	$(MAKE) -C circ constraint_bug

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
# 5x5 Square cavity for NS
NSO = $(PETSCFEM_DIR)/../petscfem-2.16g/applications/ns/ns.bin
sqcav5_dat = sqcav5.nod.tmp sqcav5.con.tmp sqcav5.fixa.tmp
sqcav5: ns sqcav5_dat
	-$(MPIRUN) -np 3 $(NSO) -case sqcav5.dat | tee sqcav5.out

sqcav5_dat: mksqcav5.m
	octave -q  mksqcav5.m

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
chunk_size_process:
	$(MAKE) -C sector chunk_size_process

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
triangle:
	$(MAKE) -C sector triangle

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
# FileStack class check
fstack:
	$(MAKE) -C fstack fstack

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
tests: runtests
	@echo "Verifying results..."
	@runtests.pl
	@echo "Done."

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
iterator: iterator.cpp
	g++ -o $@ $<

include makefile.d
