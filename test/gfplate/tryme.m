## $Id: tryme.m,v 1.1 2005/01/22 12:01:51 mstorti Exp $
source("data.m.tmp");

pref = Rgas*Tref*rhoref;
cref = sqrt(gamma*Tref*Rgas);
uini = Machin*cref;

A = [uini rhoref 0;
     0 uini 1/rhoref;
     0 rhoref*cref^2 uini];

[v,d] = eig(A);
