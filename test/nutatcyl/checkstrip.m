## Copyright (C) 2002 Mario A. Storti
##
## This file is part of Octave.
##__INSERT_LICENSE__
## $Id: checkstrip.m,v 1.1 2002/08/19 16:49:11 mstorti Exp $

## Author: Mario Storti
## Keywords: petscfem-test, viscous-force-intgrator

source("data.m.tmp");

f=aload("strip.force.tmp");
Fx = -f(1);

L = 1;
rho = 1;			# density
a = 1;				# acceleration

Vol = L*h^2;
Fx_anal = Vol * rho * a;

tol = 1e-8;			# This case is solved exactly
erro = abs(Fx-Fx_anal);
relerr = erro/Fx_anal;

printf("Fx [constant acceleration] OK ? %d\n",relerr<tol);
printf("Fx %f, Fx_anal %f, error %f, rel.error %f\n",
       Fx, Fx_anal, erro,relerr);
