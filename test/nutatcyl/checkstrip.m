## Copyright (C) 2002 Mario A. Storti
##
## This file is part of Octave.
##__INSERT_LICENSE__
## $Id: checkstrip.m,v 1.2 2003/01/10 19:17:09 mstorti Exp $

## Author: Mario Storti
## Keywords: petscfem-test, viscous-force-intgrator

source("data.m.tmp");

if !use_prisms==1
  f=aload("strip.force.tmp");
else
  f=aload("stripp.force.tmp");
endif

L = 1;
rho = 1;			# density
a = 1;				# acceleration

Vol = L*h^2;
Fx_anal = Vol * rho * a;

f_anal = [-Fx_anal, 0, 0, 0, Fx_anal*h/2, -Fx_anal*L/2];

tol = 1e-8;			# This case is solved exactly
erro = merr(f-f_anal);
relerr = erro/merr(f_anal);

printf("Fx [constant acceleration] OK ? %d\n",relerr<tol);
printf("Force-moment num. %f %f %f %f %f %f\n",f);
printf("Force-moment anal. %f %f %f %f %f %f\n",f_anal);
printf("error %f, rel.error %f\n",erro,relerr);
