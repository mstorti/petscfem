## Copyright (C) 2002 Mario A. Storti
##
## This file is part of Octave.
##__INSERT_LICENSE__
## $Id: checkstrip2d.m,v 1.2 2003/01/10 16:28:52 mstorti Exp $

## Author: Mario Storti
## Keywords: petscfem-test, viscous-force-intgrator

source("data.m.tmp");
s = 2*use_exterior_normal-1;

if use_exterior_normal
  f = aload("strip2d.force.tmp");
else 
  f = aload("strip2di.force.tmp");
endif

Fx = -f(1);

L = 1;
rho = 1;			# density
a = 1;				# acceleration

Vol = L*h;
Fx_anal = s * Vol * rho * a;

tol = 1e-8;			# This case is solved exactly
erro = abs(Fx-Fx_anal);
relerr = erro/Fx_anal;

printf("Fx [constant acceleration] OK ? %d\n",relerr<tol);
printf("Fx %f, Fx_anal %f, error %f, rel.error %f\n",
       Fx, Fx_anal, erro,relerr);
