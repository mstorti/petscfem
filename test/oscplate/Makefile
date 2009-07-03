
.PHONY: force all run tests default strip gasflow nozzle  vtubei

default: oscplate4_30deg

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
#w All oscillating plate tests (why not oscplate2b?)
tests: oscplate3 oscplate2 oscplate2b oscplate4 \
	oscplate4cos oscplate4_30deg oscplate

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
#w Oscillating plate with piecewise constant function
oscplate2: ns
	-rm -f tmp_file*.tmp oscsome2.sal
	make-depl oscplate.epl oscplate2.depl -d piecewise_fun=1
	-$(MPIEXEC) -n 1 $(NS) -case oscplate2.depl $(call in_backg,$@)
	$(OCTAVE) oscplate2.m &> oscplate2.verif.tmp

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
#w Oscillating plate with sinusoidal oscillation
oscplate: ns
	-rm -f tmp_file*.tmp oscplate1.state.tmp
	make-depl oscplate.epl oscplate.depl -d Nperiods=20 
	-$(MPIEXEC) -n 1 $(NS) -case oscplate.depl $(call in_backg,$@)
	mv -f oscplate.state.tmp oscplate1.state.tmp
	$(OCTAVE) verif1.m > oscplate1.verif.tmp

SRCS = 
PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

LOC_PROG := ns
CASE_NAME := oscplate3
PROG := $(NS_G)

STATE := oscplate.state.tmp
#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
#w Oscillating plate with sinusoidal oscillation
#w Verifies quadratic convergence of the
#w Crank Nicholson method. Runs a case with N=16,
#w 32 and 128 time steps. The solution for N=128 is
#w considered the exact, and it is checked in Octave that
#w norm(u32-u128)/norm(u16-u128) < 0.30 (it should be <0.25). 
oscplate3: ns
	-rm -f tmp_file*.tmp oscsome3.sal oscsome.sal
#
	rm -f 
	make-depl oscplate.epl oscplate3.depl -d alpha=0.5 -d N=16 
	-$(MPIEXEC) -n 1 $(NS) -case oscplate3.depl $(call in_backg,$@-1)
	mv oscplate.state.tmp oscplate_16.sal
#
	make-depl oscplate.epl oscplate3.depl -d alpha=0.5 -d N=32 
	-$(MPIEXEC) -n 1 $(NS) -case oscplate3.depl  $(call in_backg,$@-2)
	mv oscplate.state.tmp oscplate_32.sal
#
	make-depl oscplate.epl oscplate3.depl -d alpha=0.5 -d N=128 
	-$(MPIEXEC) -n 1 $(NS) -case oscplate3.depl  $(call in_backg,$@-3)
	mv oscplate.state.tmp oscplate_128.sal
#
	$(OCTAVE) oscplate3.m &> oscplate3o.sal

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
#w Oscillating plate with piecewise constant function
oscplate2b: ns
	-rm -f tmp_file*.tmp oscsome2b.sal oscsome.sal
	make-depl oscplate.epl oscplate2b.depl -d spline_fun=1  
	-$(MPIEXEC) -n 1 $(NS) -case oscplate2b.depl $(call in_backg,$@)

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
#w Oscillating plate with piecewise constant periodic function
#w `sin(t/T)' case
oscplate4: ns
	-rm -f tmp_file*.tmp oscsome4.sal oscsome.sal
	make-depl oscplate.epl oscplate4.depl -d spline_periodic_fun=1 
	-$(MPIEXEC) -n 1 $(NS) -case oscplate4.depl $(call in_backg,$@)
	mv oscsome.sal oscsome4.sal
	$(OCTAVE) oscplate4.m &> oscplate4o.sal

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
#w `cos(2*pi*t/T)' case
oscplate4cos: ns
	-rm -f tmp_file*.tmp oscsome4c.sal oscsome.sal
	make-depl oscplate.epl oscplate4c.depl \
		-d spline_periodic_fun=1  -d cos_case=1
	-$(MPIEXEC) -n 1 $(NS) -case oscplate4c.depl $(call in_backg,$@)
	mv oscsome.sal oscsome4c.sal
	$(OCTAVE) oscplate4c.m &> oscplate4co.sal

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
#w `sin(2*pi*t/T-pi/6)' case
oscplate4_30deg: ns
	-rm -f tmp_file*.tmp oscsome4_30deg.sal oscsome.sal
	make-depl oscplate.epl oscplate4_30deg.depl \
		 -d spline_periodic_fun=1 -d cos_case=0 -d sin30deg_case=1
	-$(MPIEXEC) -n 1 $(NS) -case oscplate4_30deg.depl $(call in_backg,$@)
	mv oscsome.sal oscsome4_30deg.sal
	$(OCTAVE) oscplate4_30deg.m &> oscplate4_30dego.sal


#s
