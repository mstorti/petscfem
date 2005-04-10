## $Id: mknozzle2.m,v 1.1 2005/04/10 09:14:04 mstorti Exp $
source("data.m.tmp");

pref = Rgas*Tref*rhoref;
cref = sqrt(gamma*Tref*Rgas);

uini = Machin*cref;
poutlet = pref;

w = zhomo([0 Lx1+Lx2 0 Ly],Nx1+Nx2+1,Ny+1,[1 0 1 1 0 yratio]);
[xnod,icone] = pfcm2fem(w);
icone = icone(:,[1 4 3 2]);

x = xnod(:,1);
y = xnod(:,2);

indx = find(x<Lx1);
c = 1-(DLy/Ly)*(1-abs((x-Lx1/2)/(Lx1/2)));
y(indx) = y(indx).*c(indx);

xnod = [x,y];

gplfem(xnod,icone);
gplot "malla.gpl"

return


Hin = Ly;
Hout = Ly/Aratio;
Hmean = (Hin+Hout)/2;
DH = (Hin-Hout)/2;
eta = (Ly-y)/Ly;
xi = 2*atanh(0.99)*(x-Lx/2)/Lnozzle;
eta = eta.*(Hmean-DH*tanh(xi))/Hin;
y = (Ly-eta*Hin);

xnod = [x,y];

asave("gfnozzle2.nod.tmp",xnod);
asave("gfnozzle2.con.tmp",icone);

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
pfconstr("gfnozzle2.constr-wall.tmp",[tmp,tmp],2:3,nor);

## [rho u v] at inlet
pffixa("gfnozzle2.fixa-in.tmp",inlet,1:3,[rhoref,uref,0]);

## `p' at outlet
pffixa("gfnozzle2.fixa-out.tmp",outlet,4,pref);

## `v=0' at axis
pffixa("gfnozzle2.fixa-slip.tmp",slip,3);

nnod = rows(xnod);
uini = Uref(ones(nnod,1),:);
asave("gfnozzle2.ini.tmp",uini);
