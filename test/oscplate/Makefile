# $Id: Makefile,v 1.1 2003/02/17 12:39:28 mstorti Exp $ 

.PHONY: force all run tests default strip gasflow nozzle  vtubei

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
#w All oscillating plate tests (why not oscplate2b?)
tests: oscplate3 oscplate2 oscplate2b oscplate4 \
	oscplate4cos oscplate4_30deg oscplate

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
#w Oscillating plate with piecewise constant function
oscplate2: ns
	-rm -f tmp_file*.tmp oscsome2.sal
	eperl -P -d piecewise_fun=1  oscplate.epl >oscplate2.depl
	-$(MPIRUN) -np 1 $(NS) -case oscplate2.depl >/dev/null
	$(OCTAVE) oscplate2.m &> oscplate2o.sal

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
#w Oscillating plate with sinusoidal oscillation
oscplate: ns
	-rm -f tmp_file*.tmp oscplate1.sal
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
	-rm -f tmp_file*.tmp oscsome3.sal oscsome.sal
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
	$(OCTAVE) oscplate3.m &> oscplate3o.sal

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
#w Oscillating plate with piecewise constant function
oscplate2b: ns
	-rm -f tmp_file*.tmp oscsome2b.sal oscsome.sal
	eperl -P -d spline_fun=1  oscplate.epl >oscplate2b.depl
	-$(MPIRUN) -np 1 $(NS) -case oscplate2b.depl >/dev/null

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
#w Oscillating plate with piecewise constant periodic function
#w `sin(t/T)' case
oscplate4: ns
	-rm -f tmp_file*.tmp oscsome4.sal oscsome.sal
	eperl -P -d spline_periodic_fun=1 oscplate.epl >oscplate4.depl
	-$(MPIRUN) -np 1 $(NS) -case oscplate4.depl >/dev/null
	mv oscsome.sal oscsome4.sal
	$(OCTAVE) oscplate4.m &> oscplate4o.sal

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
#w `cos(2*pi*t/T)' case
oscplate4cos: ns
	-rm -f tmp_file*.tmp oscsome4c.sal oscsome.sal
	eperl -P -d spline_periodic_fun=1  -d cos_case=1 \
			oscplate.epl >oscplate4c.depl
	-$(MPIRUN) -np 1 $(NS) -case oscplate4c.depl >/dev/null
	mv oscsome.sal oscsome4c.sal
	$(OCTAVE) oscplate4c.m &> oscplate4co.sal

#----<*>----------<*>----------<*>----------<*>----------<*>----------<*>-----
#w `sin(2*pi*t/T-pi/6)' case
oscplate4_30deg: ns
	-rm -f tmp_file*.tmp oscsome4_30deg.sal oscsome.sal
	eperl -P -d spline_periodic_fun=1 -d cos_case=0 -d sin30deg_case=1 \
			oscplate.epl >oscplate4_30deg.depl
	-$(MPIRUN) -np 1 $(NS) -case oscplate4_30deg.depl >/dev/null
	mv oscsome.sal oscsome4_30deg.sal
	$(OCTAVE) oscplate4_30deg.m &> oscplate4_30dego.sal

SRCS = 
PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

#s
