## Copyright (C) 2002 Mario A. Storti
##
## This file is part of Octave.
##__INSERT_LICENSE__
## $Id: checkrigid.m,v 1.2 2002/08/19 16:05:10 mstorti Exp $

## Author: Mario Storti
## Keywords: petscfem-test, viscous-force-intgrator

source("data.m.tmp");

f0=aload("cylinder.force.rigid.tmp");
Mz = -f0(6);

Vol = pi*RR^2*Lz;
I = Vol*RR/2;
Mz_anal = I * omega;

tol = 0.15;
erro = abs(Mz-Mz_anal);
relerr = erro/Mz_anal;

printf("Mz(rigid) OK ? %d\n",relerr<tol);
printf("Mz %f, Mz_anal %f, error %f, rel.error %f\n",
       Mz, Mz_anal, erro,relerr);
