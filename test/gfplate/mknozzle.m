## $Id: mknozzle.m,v 1.2 2005/01/31 02:35:52 mstorti Exp $
source("data.m.tmp");

pref = Rgas*Tref*rhoref;
cref = sqrt(gamma*Tref*Rgas);

uini = Machin*cref;
poutlet = pref;

w = zhomo([0 Lx 0 Ly],Nx+1,Ny+1,[1 0 1 1 0 yratio]);
[xnod,icone] = pfcm2fem(w);
icone = icone(:,[1 4 3 2]);

x = xnod(:,1);
y = xnod(:,2);

Hin = Ly;
Hout = Ly/Aratio;
Hmean = (Hin+Hout)/2;
DH = (Hin-Hout)/2;
eta = (Ly-y)/Ly;
xi = 2*atanh(0.99)*(x-Lx/2)/Lnozzle;
eta = eta.*(Hmean-DH*tanh(xi))/Hin;
y = (Ly-eta*Hin);

xnod = [x,y];

asave("gfnozzle.nod.tmp",xnod);
asave("gfnozzle.con.tmp",icone);

inlet = (1:Ny+1)';
slip = (Ny+1)*(1:Nx+1);
outlet = Nx*(Ny+1)+inlet;
wall = 1+(Ny+1)*(0:Nx)';

## Tangents at the wall
nwall = length(wall);
tan = zeros(nwall,2);
tan(2:nwall-1,:) = xnod(wall(3:nwall),:) \
    - xnod(wall(1:nwall-2),:);
tan(1,:) = [-3,4,-1]*xnod(wall(1:3),:);
tan(nwall,:) = [1,-4,3]*xnod(wall(nwall+(-2:0)),:);

## Normals
nor = [-tan(:,2),tan(:,1)];

xwall = xnod(wall,:);
xwall2 = xwall+0.1*nor;

## Slip at the wall
tmp = complement(inlet,wall)';
pfconstr("gfnozzle.constr-wall.tmp",[tmp,tmp],2:3,nor);

## [rho u v] at inlet
pffixa("gfnozzle.fixa-in.tmp",inlet,1:3,[rhoref,uref,0]);

## `p' at outlet
pffixa("gfnozzle.fixa-out.tmp",outlet,4,pref);

## `v=0' at axis
pffixa("gfnozzle.fixa-slip.tmp",slip,3);

nnod = rows(xnod);
uini = Uref(ones(nnod,1),:);
asave("gfnozzle.ini.tmp",uini);
