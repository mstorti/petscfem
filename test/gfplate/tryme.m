## $Id: tryme.m,v 1.2 2005/01/22 22:10:20 mstorti Exp $
source("data.m.tmp");

pref = Rgas*Tref*rhoref;
cref = sqrt(gamma*Tref*Rgas);
uini = Machin*cref;

A = [uini rhoref 0;
     0 uini 1/rhoref;
     0 rhoref*cref^2 uini];

[v,la] = eig(A);
la = diag(la);
