# $Id: Makefile,v 1.81 2004/11/06 13:22:27 mstorti Exp $ 

.PHONY: all run  compile lclean laplace adv nso laplace_tests force	\
        ctest gprof dxhook dx_tests dx_hook_test_post_run dx		\
	dx_hook_test_post_run_file_record test_read_ini 		\
	test_disconnected_graph

BOPT := g_c++
default: 
	$(MAKE) case=iisd_sbp1 np=2 runs

isp_test:
	$(MAKE) CASE=qharm np=1 case=iisd_sbp2 run

###--------<*>--------<*>--------<*>--------<*>--------<*>--------<*>--------<*>
tests: test_re1000 test_iisd test_iisd2 test_meas_perf		\
	test_g_body test_fs test_zero_w_proc test_zero_w_proc2	\
	test_triangle dx_tests test_qharm_3d test_shell_hook 	\
	test_read_ini test_disconnected_graph

test_shell_hook: test_shell_hook1 test_shell_hook2

test_shell_hook1:
	$(MAKE) CASE=qharm case=shell_hook run
	-rm -f outvector0.out

test_shell_hook2:
	$(MAKE) CASE=qharm case=shell_hook2 run
	-rm -f outvector0.out

#---<*>---//---<*>---//---<*>---//---<*>---//---<*>---// 
#w Tests the adaptor_pg class when ndimel<ndim (quad panels in 3D). 
test_qharm_3d:
# A Laplace problem in surfaces is considered. We simply consider a 2D
# plane but simulating it as a plane in space. First, for reference,
# we solve the strictly 2D problem, and then we solve for the plane
# passing through the `y' axis at several angles with respect to the
# `xy' plane. (0,45 and 90deg).
	-rm -f qharm.state_g3d_ref.tmp
	$(MAKE) CASE=qharm case=g3d_ref run
	-cp outvector0.out qharm.state_g3d_ref.tmp
	-rm -f qharm.state_g3d0.tmp
	$(MAKE) CASE=qharm case=g3d_0 run
	-cp outvector0.out qharm.state_g3d_0.tmp
	-rm -f qharm.state_g3d45.tmp
	$(MAKE) CASE=qharm case=g3d_45 run
	-cp outvector0.out qharm.state_g3d_45.tmp
	-rm -f qharm.state_g3d90.tmp
	$(MAKE) CASE=qharm case=g3d_90 run
	-cp outvector0.out qharm.state_g3d_90.tmp
	octave -qH check_g3d.m > verif.g3d.tmp

dx_tests_np1 := dx_hook_test_post_run_file_record	\
		dx_hook_test_post_run_files		\
		dx_hook_test_synchro			\
		dx_hook_test_allf
dx_tests_np2 := $(patsubst dx_hook_test_%,dx_test_%_np2,$(dx_tests_np1))
dx_tests_np1: $(dx_tests_np1) 
dx_tests_np2: $(dx_tests_np2) 
dx_tests: dx_tests_np1 dx_tests_np2

dx_test_%_np2: 
	$(MAKE) np=2 dx_hook_test_$*

dx:
	$(MAKE) -C ../../dx

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

LOC_PROG := ns
CASE_NAME := sqcav
PROG := $(NS)

###--------<*>--------<*>--------<*>--------<*>--------<*>--------<*>--------<*> 
NP = 1

MCH = -machinefile ./machi.dat

NP := $(shell procsel $(PETSCFEM_PROCTABLE))
ifeq ($(np),)
np := 1
endif
output_file := output.CASE_$(CASE).np_$(np).case_$(case)$(meas_perf).out.tmp

#w Runs program
run:  data
	-${MPIRUN} -np $(np) $(MCH) $(PROG) -case $(CASE).depl \
	 	      	  > $(output_file)

export CASE
runs: CASE := sqcav
runs:
	$(MAKE) run

runs2: CASE := sqcav2
runs2:
	$(MAKE) run

runq: CASE := qharm
runq: force
#	$(MAKE) case=iisd_sbp1 run
	for n in 1 2 4 6 8 16 32 64 128 ; \
		do time $(MAKE) case=iisd_sbp$$n run ; \
		mv output.out output.iisd$$n.out ; \
	done

runq1: CASE := qharm
runq1: 
	$(MAKE) case=iisd_sbp1 run

runs: CASE := sqcav
test_big:
	time $(MAKE) case=big runs

test_g_body:
	$(MAKE) case=g_body runs
	$(OCTAVE) < check.g_body.m > check.g_body.verif.tmp

###--------<*>--------<*>--------<*>--------<*>--------<*>--------<*>--------<*> 
# tests the square cavity in a 20x20 mesh at Re=1000
# with respect with a reference solution obtained also with Petscfem and
# verified with Ghia results
#w IISD for square cavity test in 20x20 mesh on one processor
test_re1000:
	$(MAKE) case=weak_form_0 runs
	$(MAKE) case=weak_form_1 runs
	$(OCTAVE) < check_re1000.m > check.re1000.verif.tmp

###--------<*>--------<*>--------<*>--------<*>--------<*>--------<*>--------<*> 
#w IISD for square cavity test in 20x20 mesh on two processors
test_iisd2: 
	$(MAKE) case=lu np=1 runs
	$(MAKE) case=iisd_sbp1 np=2 runs
	$(MAKE) case=iisd_sbp4 np=2 runs
	$(MAKE) case=iisd_sbp16 np=2 runs
	$(OCTAVE) < check_iisd2.m > check.iisd.verif.np2.tmp

REDIR := >
# REDIR := | tee

###--------<*>--------<*>--------<*>--------<*>--------<*>--------<*>--------<*> 
#w Tests processors with zero eqs./elements (zero weight)
test_zero_w_proc: CASE := sqcav
test_zero_w_proc:
# 	Runs with one and two process and compares
#
#	One process:
	$(MAKE) case=zwproc_ref np=1 data
	echo $(HOSTNAME) > machi.dat
	echo 1. > weights.dat
	-${MPIRUN} -np 1 -machinefile ./machi.dat $(PROG) -case $(CASE).depl \
	$(REDIR) output.zwproc_ref.out.tmp
#	Two process:
	$(MAKE) case=zwproc np=2 data
	echo $(HOSTNAME) > machi.dat
	echo 0. > weights.dat
	echo 1. >> weights.dat
	-${MPIRUN} -np 2 -machinefile ./machi.dat $(PROG) -case $(CASE).depl \
	$(REDIR) output.zwproc.out.tmp
#
	$(OCTAVE) check_zwproc.m > check.zwproc.tmp

###--------<*>--------<*>--------<*>--------<*>--------<*>--------<*>--------<*> 
#w Tests processors with zero eqs./elements (zero weight)
test_zero_w_proc2: CASE := sqcav
test_zero_w_proc2:
# 	Runs with one and two process and compares
#
#	One process:
	$(MAKE) case=zwproc2_ref np=1 data
	echo $(HOSTNAME) > machi.dat
	echo 1. >> weights.dat
	-${MPIRUN} -np 1 -machinefile ./machi.dat $(PROG) -case $(CASE).depl \
	$(REDIR) output.zwproc2_ref.out.tmp
#
#	Two process:
	$(MAKE) case=zwproc2 np=2 data
	echo $(HOSTNAME) > machi.dat
	echo 0. > weights.dat
	echo 1. >> weights.dat
	-${MPIRUN} -np 2 -machinefile ./machi.dat $(PROG) -case $(CASE).depl \
	$(REDIR) output.zwproc2.out.tmp
#
	$(OCTAVE) check_zwproc2.m > check.zwproc2.tmp

test_meas_perf:
	$(MAKE) case=lu np=1 meas_perf=1 runs

###--------<*>--------<*>--------<*>--------<*>--------<*>--------<*>--------<*>
#w Tests IISD solver on the square cavity problem
test_iisd:
	$(MAKE) case=lu np=1 runs
	$(MAKE) case=iisd_sbp1 runs
	$(MAKE) case=iisd_sbp2 runs
	$(MAKE) case=iisd_sbp4 runs
	$(MAKE) case=iisd_sbp8 runs
	$(MAKE) case=iisd_sbp16 runs
	$(OCTAVE) < check_iisd.m > check.iisd.verif.np1.tmp

###--------<*>--------<*>--------<*>--------<*>--------<*>--------<*>--------<*>
#w Tests square cavity at Re=400 with fractional step solver. 
test_fs:
	$(MAKE) case=fractional_step_re400 runs2
	$(OCTAVE) < check_fs.m > check.fs.verif.tmp

###--------<*>--------<*>--------<*>--------<*>--------<*>--------<*>--------<*>
#w Tests triangles
test_triangle:
	$(MAKE) case=weak_form_1 use_triangles=1 runs

###--------<*>--------<*>--------<*>--------<*>--------<*>--------<*>--------<*> 
output_file1 := output.CASE_sqcav.np_$(np).case_dx_file1.out.tmp
#w Checks reading states from file in severa records. 
dx_hook_test_post_run_file_record: dx
#
# A series of 20 states are generated in file 'sqcav.state_0.tmp'. 
# Then simulates that dx is running and reads 11 states from differents
# records in the same file. Then stops by sending record=-1 to DX. 
#
# Clean files
	-rm -f dx_output.file1.tmp $(output_file1)
# Generates all states in a file 'sqcav.state_0.tmp'
	$(MAKE) case=dx_file0 runs
# Runs DX in background and PETSc-FEM 
	dx -mdf ./epimport.mdf -cache off -memory 100 \
			-script qq1.net > dx_output.file1.tmp &
	$(MAKE) case=dx_file1 runs
	echo "------<snip>-----------" >> $(output_file1)
	cat dx_output.file1.tmp >> $(output_file1)

###--------<*>--------<*>--------<*>--------<*>--------<*>--------<*>--------<*> 
output_file3 := output.CASE_sqcav.np_$(np).case_dx_file3.out.tmp
#w Checks reading states from several files
dx_hook_test_post_run_files: dx
#
# A series of 20 states are generated in files 'sqcav.state_<n>.tmp'. 
# Then simulates that dx is running and reads several states from different
# files. Then stops by sending record=-1 to DX. 
#
# Clean files
	-rm -f dx_output.file3.tmp $(output_file3) sqcav.state_*.tmp
# Generates all states in a file 'sqcav.state_0.tmp'
	$(MAKE) case=dx_file2 runs
# Runs DX in background and PETSc-FEM 
	dx -mdf ./epimport.mdf -cache off -memory 100 \
			-script qq3.net > dx_output.file3.tmp &
	$(MAKE) case=dx_file3 runs
	echo "------<snip>-----------" >> $(output_file3)
	cat dx_output.file3.tmp >> $(output_file3)

###--------<*>--------<*>--------<*>--------<*>--------<*>--------<*>--------<*> 
output_synchro := output.CASE_sqcav.np_$(np).case_dx_synchro.out.tmp
dx_output := dx_output.synchro.tmp
#w Checks reading states in synchronous form
dx_hook_test_synchro: dx
#
# DX and PETSc-FEM run at the same time, exchanging information
#
# Clean files
	-rm -f $(dx_output) $(output_synchro)
# Runs DX in background and PETSc-FEM 
	dx -mdf ./epimport.mdf -cache off -memory 100 \
			-script qq_synchro.net > $(dx_output) &
	$(MAKE) case=dx_synchro runs
	echo "------<snip>-----------" >> $(output_synchro)
	cat $(dx_output) >> $(output_synchro)

###--------<*>--------<*>--------<*>--------<*>--------<*>--------<*>--------<*> 
output_allf := output.CASE_sqcav.np_$(np).case_dx_allf.out.tmp
dx_output := dx_output.allf.tmp
#w Checks reading states in synchronous form
dx_hook_test_allf: dx
#
# DX and PETSc-FEM run at the same time, exchanging information
#
# Clean files
	-rm -f $(dx_output) $(output_allf)
# Runs DX in background and PETSc-FEM 
	dx -mdf ./epimport.mdf -cache off -memory 100 \
			-script qq_synchro.net > $(dx_output) &
	$(MAKE) case=dx_allf runs
	echo "------<snip>-----------" >> $(output_allf)
	cat $(dx_output) >> $(output_allf)

data: force
	$(MAKE) $(CASE).depl
	if [ -f mk$(CASE).m ] ; then $(OCTAVE) mk$(CASE).m ; fi

# EPERLFLAGS = -d case=$$case
$(CASE).depl: force
#	echo in target \$(CASE).depl : CASE $(CASE)

local_clean::
	-chmod +w *.depl
	-rm *.depl *.sal tmp_file*.tmp &> /dev/null

shell_hook_%:
	@echo target $@, step: $$petscfem_step, time: $$petscfem_time

###--------<*>--------<*>--------<*>--------<*>--------<*>--------<*>--------<*> 
test_read_ini:
	$(MAKE) case=read_ini np=4 runs
	octave -qH check.read_ini.m > verif.read_ini.tmp

###--------<*>--------<*>--------<*>--------<*>--------<*>--------<*>--------<*> 
test_press_term:
	$(MAKE) case=press_term runs

###--------<*>--------<*>--------<*>--------<*>--------<*>--------<*>
#w Checks for disconnected graphs. 
test_disconnected_graph:
	$(MAKE) case=disconnected_graph CASE=simple run

###--------<*>--------<*>--------<*>--------<*>--------<*>--------<*>
#w Makes local link for running with debugger
link:
	$(MAKE) -C $(PETSCFEM_DIR) $(LOC_PROG)
	ln -f $(PROG) $(LOC_PROG)
