## $Id: mkgfshock.m,v 1.2 2005/02/25 21:25:33 mstorti Exp $
source("data.m.tmp");

w = zhomo([0 Lx/Nx 0 Lx ],2,Nx+1);
[xnod,icone] = pfcm2fem(w);
xnod = xnod(:,[2 1]);
nnod = size(xnod,1);
x = xnod(1:nnod,1);

asave("gfshock.nod.tmp",xnod);
asave("gfshock.con.tmp",icone);

## rho,u,v at inlet
inlet = [1];
pffixa("gfshock.fixa-in.tmp", \
       inlet,1:3,[rhoin uin 0]);

## p at outlet
outlet = [Nx+1];
pffixa("gfshock.fixa-outlet.tmp",outlet,4,pout)

## periodic
pfperi("gfshock.peri.tmp",Nx+1+(1:Nx+1)',(1:Nx+1)',1:4);

## Fixa on reference nodes
Uini = [rhoout,uout,0,pout];
Uini = Uini(ones(nnod,1),:);

asave("gfshock.ini.tmp",Uini);
