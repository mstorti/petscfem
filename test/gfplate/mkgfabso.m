## $Id: mkgfabso.m,v 1.16 2005/01/23 22:46:34 mstorti Exp $
source("data.m.tmp");

poutlet = pref;

dx = Lx/Nx;
w = zhomo([0 dx 0 Lx ],2,Nx+1);
[x2,ico2] = pfcm2fem(w);
x2 = x2(:,[2 1]);
nnod = size(x2,1);

[xnod,icone] = extrude(x2,ico2,1,dx);
x = xnod(:,1);

Oy = eye(3);
Oy([3,1],[3,1]) = [cos(theta),sin(theta);
		   -sin(theta),cos(theta)];
Oz = eye(3);
Oz(1:2,1:2) = [cos(phi),-sin(phi);
	       sin(phi),cos(phi)];

Orot = Oz*Oy; # Rotation matrix
xnod = xnod*Orot';

ndim=3;
ndof=ndim+2;

## Periodic b.c.'s
fid = fopen("gfabso.peri.tmp","w");
for k=1:Nx+1
  for dof=1:ndof
    for l=1:3
      node = k+l*(Nx+1);
      fprintf(fid,"%g %d %d    %g %d %d\n",
	      -1,node,dof,1,k,dof);
    endfor
  endfor
endfor
fclose(fid);

nnod = rows(xnod);
nnod==4*(Nx+1) || error("inconsistency...");
## Fictitious nodes at outlet 
## ... 4*Nx+4 nnod+1 nnod+2
## Fictitious nodes at inlet
## nnod+3, nnod+4, 1
xnod = [xnod;
	xnod(nnod,:);
	xnod(nnod,:);
	xnod(1,:);
	xnod(1,:)];

asave("gfabso.nod.tmp",xnod);
asave("gfabso.con.tmp",icone);

## Absorbing b.c.'s
abso1 = [nnod+(0:-1:-2), nnod+[1,2]];
asave("gfabso.con-abso1.tmp",abso1);

abso0 = [1:3,nnod+[3,4]];
asave("gfabso.con-abso0.tmp",abso0);

## Fixa on reference nodes
Uref = [rhoref,[uref,0,0]*Orot',pref];
ref = [nnod+2,nnod+4];
pffixa("gfabso.fixa-ref.tmp",ref,1:5,Uref)

asave("gfabso.some-nodes.tmp",(1:Nx+1)');

nnod2 = size(xnod,1);
Uini = Uref(ones(nnod2,1),:);
Uini(nnod+[1,3],:) = 0;		# lagrange multipliers to 0

dw = 0.2*[0 0 1 0 0]; ## perturbation
dw(2:4) = dw(2:4)*Orot';

dfx = exp(-((x-Lx/2)/sigma).^2);
Uini(1:nnod,:) = Uini(1:nnod,:) + dfx*dw;

asave("gfabso.ini.tmp",Uini);
