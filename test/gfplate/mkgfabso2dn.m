## $Id: mkgfabso2dn.m,v 1.13 2005/02/04 12:22:21 mstorti Exp $
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
pffixa("gfabso2dn.fixa-in.tmp", \
       inlet,1:3,[rhoref [uref 0]*Orot]);

## p at outlet
outlet = [Nx+1;2*Nx+2];
pffixa("gfabso2dn.fixa-outlet.tmp", \
       outlet,4,pref);

## other
t = [0,1]*Orot;
slip = (1:nnod)';
pfconstr("gfabso2dn.fixa-slip.tmp",[slip,slip],2:3,t);

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

xe = pfnd2ele(xnod,icone,xnod(:,1));
Gb = zeros(Nx,2);
lgb = Lx/4; 			# Length where Gb is applied
xi = (xe-Lx/2)/(lgb/2);
Gb(:,1) = (1-xi.^2).*(abs(xi)<1);

fid = fopen("gfabso2dn.con.tmp","w");
for k=1:Nx
  fprintf(fid,"%d %d %d %d    %g %g\n",
	  icone(k,:),Gb(k,:));
endfor
fclose(fid);

asave("gfabso2dn.nod.tmp",xnod);
## asave("gfabso2dn.con.tmp",icone);

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

nor = [1,0]*Orot;
fid = fopen("gfabso2dn.con-nabso.tmp","w");
fprintf(fid,"%d %d %d     %g %g\n",1,nnod+[5,6],-nor);
fprintf(fid,"%d %d %d     %g %g\n",Nx+2,nnod+[7,8],-nor);
fprintf(fid,"%d %d %d     %g %g\n",Nx+1,nnod+[1,2],nor);
fprintf(fid,"%d %d %d     %g %g\n",2*Nx+2,nnod+[3,4],nor);
fclose(fid);

## Fixa on reference nodes
Uref = [rhoref,[uref,0]*Orot,pref];
ref = [Nx+1;2*Nx+2];
pffixa("gfabso2dn.fixa-ref.tmp",nnod+[2,4,6,8],1:4,Uref);

## Fix all fields at outlet and fictitious values
## nodes = [1,Nx+1,Nx+2,2*Nx+2,nnod+(1:8)]';
twall_con = [1,nnod+5;
	     Nx+2,nnod+7;
	     Nx+1,nnod+1;
	     2*Nx+2,nnod+3];
asave("gfabso2dn.twall-con.tmp",twall_con);
pffixa("gfabso2dn.fixa-twall.tmp", \
       nnod+(1:2:7)',2,Tref*[0.9,0.9,1.1,1.1]')
pffixa("gfabso2dn.fixa-unused.tmp", \
       nnod+(1:2:7)',[3,4]);
pffixa("gfabso2dn.fixa-u.tmp",
       [1,Nx+2,Nx+1,2*Nx+2]',[2,3]);
asave("gfabso2dn.some-nodes.tmp",(1:nnod/2)');

nnod2 = size(xnod,1);
Uini = Uref(ones(nnod2,1),:);
Uini(nnod+[1:2:7],:) = 0;	# lagrange multipliers to 0

du=0.0;
dw = du*[0 0 0 1]; ## perturbation for all waves
## dw = du*[1 0 0 0]; ## entropy wave
## dw = du*[1 cref/rhoref 0 cref^2]; ## forward acoustic wave
## dw = du*[-1 cref/rhoref 0 -cref^2]; ## backward acoustic wave
dw(2:3) = dw(2:3)*Orot;

dfx = exp(-((x-Lx/2)/sigma).^2);
Uini(1:nnod,:) = Uini(1:nnod,:) + dfx*dw;

asave("gfabso2dn.ini.tmp",Uini);
