## $Id: mkgfabso.m,v 1.3 2005/01/20 22:50:35 mstorti Exp $
source("data.m.tmp");

pref = Rgas*Tref*rhoref;
cref = sqrt(gamma*Tref*Rgas);

uini = Machin*cref;
poutlet = pref;

w = zhomo([0 Lx/Nx 0 Lx ],2,Nx+1);
[xnod,icone] = pfcm2fem(w);
xnod = xnod(:,[2 1]);

## rho,u,v at inlet
inlet = [1;Nx+2];
pffixa("gfabso.fixa-in.tmp",inlet,1:3,[rhoref uini 0.0])

## rho,u,v at inlet
outlet = [Nx+1;2*Nx+2];
pffixa("gfabso.fixa-outlet.tmp",outlet,4,pref)

## other
slip = (1:2*Nx+2)';
slip = complement(inlet,slip);
slip = complement(outlet,slip);
pffixa("gfabso.fixa-slip.tmp",outlet,3)

asave("gfabso.nod.tmp",xnod);
asave("gfabso.con.tmp",icone);

indx = (Nx-3:Nx+1);
asave("gfabso.con-abso1.tmp",
      [indx;
       Nx+1+indx]);

asave("gfabso.some-nodes.tmp",(1:Nx+1)');

Uini = [rhoref,uini,0,pref];
nnod = size(xnod,1);
Uini = Uini(ones(nnod,1),:);

x = xnod(:,1);
Uini(:,2) = Uini(:,2) + du*exp(-((x-Lx/2)/sigma).^2);

asave("gfabso.ini.tmp",Uini);

