## $Id: mkgfperi.m,v 1.1 2005/02/14 04:31:31 mstorti Exp $
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

## Periodic in `y' direction
slip = (1:nnod)';
y0 = (1:(Nx+1))';
y1 = Nx+1+(1:(Nx+1))';
pfperi("gfperi.peri-y.tmp",y1,y0,1:4);

## Periodic in `x' direction
pfperi("gfperi.peri-x.tmp",Nx+1,1,1:4);

asave("gfperi.nod.tmp",xnod);
asave("gfperi.con.tmp",icone);

nnod2 = size(xnod,1);
Uini = Uref(ones(nnod2,1),:);
Uini(nnod+[1:2:7],:) = 0;	# lagrange multipliers to 0

du=0.3;
dw = du*[0 0 0 1]; ## perturbation for all waves
## dw = du*[1 0 0 0]; ## entropy wave
## dw = du*[1 cref/rhoref 0 cref^2]; ## forward acoustic wave
## dw = du*[-1 cref/rhoref 0 -cref^2]; ## backward acoustic wave
dw(2:3) = dw(2:3)*Orot;

dfx = exp(-((x-Lx/2)/sigma).^2);

Uini = U1(ones(nnod2,1),:);
dw = U2-U1;
xx = xnod(:,1);
dfx = (xx>0.3*L & xx<0.6*L);
Uini = Uini + dfx*dw;

Uini(1:nnod,:) = Uini(1:nnod,:) + dfx*dw;

asave("gfperi.ini.tmp",Uini);
asave("gfperi.some-nodes.tmp",(1:nnod/2)');
