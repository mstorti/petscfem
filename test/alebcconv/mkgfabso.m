## $Id: mkgfabso.m,v 1.10 2006/06/21 11:33:08 mstorti Exp $
source("data.m.tmp");

dx = Lx/Nx;
w = zhomo([0 dx 0 Lx ],2,Nx+1);
[xnod,icone] = pfcm2fem(w);
xnod=xnod(:,[2,1]);

pfperi("gfabso.peri.tmp",
       (1:Nx+1)',Nx+1+(1:Nx+1)',(1:4));
nnod = rows(xnod);

asave("gfabso.nod.tmp",[xnod,xnod]);
asave("gfabso.con.tmp",icone);

pffixa("gfabso.fixa.tmp",[1;Nx+1],(2:3),[vmesh,0]);

uini = [rhoref,ufl,0,pref];
uini = uini(ones(nnod,1),:);
asave("gfabso.ini.tmp",uini);

fixed = [1,3,4];
pffixa("gfabso.fixa-dofs.tmp",(1:Nx+1)',fixed,Uref(fixed));
