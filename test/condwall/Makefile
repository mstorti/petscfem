
default: condwall

# BOPT = O_c++
PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

LOC_PROG := ns
PROG := $(NS)
CASE := condwall

MCH = -nolocal -machinefile ./machi.dat
NP := $(shell procsel -n ./proctable)

.PHONY: default condwall

condwall:
	$(MAKE) run

###--------<*>--------<*>--------<*>--------<*>--------<*>--------<*>
#w Makes local link for running with debugger
movwall.efn: movwall.cpp $(LIBPETSCFEM)

IN_BACKG := </dev/null >burgers.$(CASE).out.tmp
# IN_BACKG := 

#w Runs program
run: movwall.efn
	mkstepsdir.pl
	make-depl $(CASE).epl $(CASE).depl
	-${MPIEXEC} $(MCH) -n $(NP) $(PROG) -case $(CASE).depl \
		$(IN_BACKG)

#w Makes some cleaning
cw_shell_hook_init:

cw_shell_hook_pre:

cw_shell_hook_post: zip

cw_shell_hook_close:

###--------<*>--------<*>--------<*>--------<*>--------<*>--------<*>
#w Runs PETSc-FEM instance for DX
dx:
	make-depl $(CASE).epl $(CASE)-dx.depl -d dx=1
	-${MPIEXEC} -n 1 $(PROG) -case $(CASE)-dx.depl

zip:
	@shop -s nullglob ;				\
	for f in STEPS/$(CASE).state_*.tmp ;		\
	do echo gzipping $$f ; gzip -f $$f ; done

dx_make_command:
	octave -qH fixstate.m

###--------<*>--------<*>--------<*>--------<*>--------<*>--------<*>
#w Makes local link for running with debugger
link: movwall.efn
	$(MAKE) -C $(PETSCFEM_DIR) $(LOC_PROG)
	ln -f $(PROG) $(LOC_PROG)

yuv_compress:
	yuv=YUV/$(CASE).$$dx_step.yuv ;		\
	tiff=$(CASE).tiff ;			\
	if [ -f $$tiff ] ; then			\
		if [ -f $$yuv ] ;		\
			then rm -f $$yuv ;	\
		fi ;				\
		convert $$tiff $$yuv ;		\
		gzip -f $$yuv ;			\
		rm -f $$tiff ;			\
	fi
