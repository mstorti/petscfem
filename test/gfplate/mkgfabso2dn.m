## $Id: mkgfabso2dn.m,v 1.3 2005/01/27 14:43:53 mstorti Exp $
source("data.m.tmp");

poutlet = pref;

w = zhomo([0 Lx/Nx 0 Lx ],2,Nx+1);
[xnod,icone] = pfcm2fem(w);
xnod = xnod(:,[2 1]);
nnod = size(xnod,1);
x = xnod(1:nnod,1);

Orot = [cos(alpha),sin(alpha);
	-sin(alpha),cos(alpha)]; # Rotation matrix
xnod = xnod*Orot;

## rho,u,v at inlet
inlet = [1;Nx+2];
pffixa("gfabso2dn.fixa-in.tmp",inlet,1:3,[rhoref [uref 0]*Orot])

## p at outlet
outlet = [Nx+1;2*Nx+2];
pffixa("gfabso2dn.fixa-outlet.tmp",outlet,4,pref)

## other
slip = (1:nnod)';
pffixa("gfabso2dn.fixa-slip.tmp",slip,1+longindx)

## Fictitious nodes at outlet 
## ... 2*Nx+2 nnod+3 nnod+4
## ... Nx+1   nnod+1 nnod+2
## nnod+5, nnod+6, 1
## nnod+7, nnod+8, Nx+2
xnod = [xnod;
	xnod(Nx+1,:);
	xnod(Nx+1,:);
	xnod(2*Nx+2,:);
	xnod(2*Nx+2,:);
	xnod(1,:);
	xnod(1,:);
	xnod(Nx+2,:);
	xnod(Nx+2,:)];

asave("gfabso2dn.nod.tmp",xnod);
asave("gfabso2dn.con.tmp",icone);

## Absorbing b.c.'s
abso1 = [Nx+1:-1:Nx-1 nnod+[1,2];
	 2*Nx+2+(0:-1:-2),nnod+[3,4]];
asave("gfabso2dn.con-abso1.tmp",abso1);

abso0 = [1:3,nnod+[5,6];
	 Nx+(2:4),nnod+[7,8]];
asave("gfabso2dn.con-abso0.tmp",abso0);

## New absorbing b.c.'s
abso1 = [Nx+1 nnod+[1,2];
	 2*Nx+2,nnod+[3,4]];
asave("gfabso2dn.con-nabso1.tmp",abso1);

abso0 = [1,nnod+[5,6];
	 Nx+2,nnod+[7,8]];
asave("gfabso2dn.con-nabso0.tmp",abso0);

## Fixa on reference nodes
Uref = [rhoref,[uref,0]*Orot,pref];
ref = [Nx+1;2*Nx+2];
pffixa("gfabso2dn.fixa-ref.tmp",nnod+[2,4,6,8],1:4,Uref)

asave("gfabso2dn.some-nodes.tmp",(1:nnod/2)');

nnod2 = size(xnod,1);
Uini = Uref(ones(nnod2,1),:);
Uini(nnod+[1:2:7],:) = 0;	# lagrange multipliers to 0

dw = 0.02*[1 0 0 0]; ## entropy wave
## dw = drho*[1 cref/rhoref 0 cref^2]; ## forward acoustic wave
## dw = drho*[-1 cref/rhoref 0 -cref^2]; ## backward acoustic wave
dw(2:3) = dw(2:3)*Orot;

dfx = exp(-((x-Lx/2)/sigma).^2);
Uini(1:nnod,:) = Uini(1:nnod,:) + dfx*dw;

asave("gfabso2dn.ini.tmp",Uini);
