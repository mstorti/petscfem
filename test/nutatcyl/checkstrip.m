## Copyright (C) 2002 Mario A. Storti
##
## This file is part of Octave.
##__INSERT_LICENSE__
## $Id: checkstrip.m,v 1.3 2003/01/11 00:18:25 mstorti Exp $

## Author: Mario Storti
## Keywords: petscfem-test, viscous-force-intgrator

source("data.m.tmp");

name = "strip";
if use_prisms, name = "stripp"; endif
if !use_exterior_normal, name = [name "i"]; endif
f=aload([name ".force.tmp"]);

L = 1;
rho = 1;			# density
a = 1;				# acceleration

Vol = L*h^2;
Fx_anal = Vol * rho * a;

f_anal = -(2*use_exterior_normal-1)*\
    [-Fx_anal, 0, 0, 0, Fx_anal*h/2, -Fx_anal*L/2];

tol = 1e-8;			# This case is solved exactly
erro = merr(f-f_anal);
relerr = erro/merr(f_anal);

printf("Fx [constant acceleration] OK ? %d\n",relerr<tol);
printf("Force-moment num. %f %f %f %f %f %f\n",f);
printf("Force-moment anal. %f %f %f %f %f %f\n",f_anal);
printf("error %f, rel.error %f\n",erro,relerr);
