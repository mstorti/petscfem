## $Id: mkgfabso2dn.m,v 1.16.2.1 2005/02/14 04:33:46 mstorti Exp $
source("data.m.tmp");

poutlet = pref;

mm = 1;
w = zhomo([0 mm*Lx/Nx 0 Lx ],mm,Nx+1);
[xnod,icone] = pfcm2fem(w);
xnod = xnod(:,[2 1]);
nnod = size(xnod,1);
x = xnod(1:nnod,1);

Orot = [cos(alpha),sin(alpha);
	-sin(alpha),cos(alpha)]; # Rotation matrix
xnod = xnod*Orot;

## rho,u,v at inlet
inlet = [1;Nx+2];
pffixa("gfmovshock.fixa-in.tmp", \
       inlet,1:3,[rhoref [uref 0]*Orot]);

## p at outlet
outlet = [Nx+1;2*Nx+2];
pffixa("gfmovshock.fixa-outlet.tmp", \
       outlet,4,pref);

## other
t = [0,1]*Orot;
slip = (1:nnod)';
pfconstr("gfmovshock.fixa-slip.tmp",[slip,slip],2:3,t);

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

fid = fopen("gfmovshock.con.tmp","w");
for k=1:Nx
  fprintf(fid,"%d %d %d %d    %g %g\n",
	  icone(k,:),Gb(k,:));
endfor
fclose(fid);

asave("gfmovshock.nod.tmp",xnod);
## asave("gfmovshock.con.tmp",icone);

## Absorbing b.c.'s
abso1 = [Nx+1:-1:Nx-1 nnod+[1,2];
	 2*Nx+2+(0:-1:-2),nnod+[3,4]];
asave("gfmovshock.con-abso1.tmp",abso1);

abso0 = [1:3,nnod+[5,6];
	 Nx+(2:4),nnod+[7,8]];
asave("gfmovshock.con-abso0.tmp",abso0);

## New absorbing b.c.'s
abso1 = [Nx+1 nnod+[1,2];
	 2*Nx+2,nnod+[3,4]];
asave("gfmovshock.con-nabso1.tmp",abso1);

abso0 = [1,nnod+[5,6];
	 Nx+2,nnod+[7,8]];
asave("gfmovshock.con-nabso0.tmp",abso0);

nor = [1,0]*Orot;
fid = fopen("gfmovshock.con-nabso.tmp","w");
fprintf(fid,"%d %d %d     %g %g\n",1,nnod+[5,6],-nor);
fprintf(fid,"%d %d %d     %g %g\n",Nx+2,nnod+[7,8],-nor);
fprintf(fid,"%d %d %d     %g %g\n",Nx+1,nnod+[1,2],nor);
fprintf(fid,"%d %d %d     %g %g\n",2*Nx+2,nnod+[3,4],nor);
fclose(fid);

## Periodic b.c.'s at ends
pfperi("gfabso2dn.peri.tmp",[Nx+2,2*Nx+2]',[1,Nx+2]',1:4);

## Fixa on reference nodes
Uref = [rhoref,[uref,0]*Orot,pref];
ref = [Nx+1;2*Nx+2];
pffixa("gfmovshock.fixa-ref.tmp",nnod+[2,4,6,8],1:4,Uref);

## Fix all fields at outlet and fictitious values
## nodes = [1,Nx+1,Nx+2,2*Nx+2,nnod+(1:8)]';
twall_con = [1,nnod+5;
	     Nx+2,nnod+7;
	     Nx+1,nnod+1;
	     2*Nx+2,nnod+3];
asave("gfmovshock.twall-con.tmp",twall_con);
pffixa("gfmovshock.fixa-twall.tmp", \
       nnod+(1:2:7)',2,Tref*[0.9,0.9,1.1,1.1]')
pffixa("gfmovshock.fixa-unused.tmp", \
       nnod+(1:2:7)',[3,4]);
pffixa("gfmovshock.fixa-u.tmp",
       [1,Nx+2,Nx+1,2*Nx+2]',[2,3]);
asave("gfmovshock.some-nodes.tmp",(1:nnod/2)');

nnod2 = size(xnod,1);
Uini = Uref(ones(nnod2,1),:);
Uini(nnod+[1:2:7],:) = 0;	# lagrange multipliers to 0

if 0
  du=0.2;
  dw = du*[0 0 0 1]; ## perturbation for all waves
  ## dw = du*[1 0 0 0]; ## entropy wave
  ## dw = du*[1 cref/rhoref 0 cref^2]; ## forward acoustic wave
  ## dw = du*[-1 cref/rhoref 0 -cref^2]; ## backward acoustic wave
  dw(2:3) = dw(2:3)*Orot;
endif

M1 = 6;
p1 = 1;
rho1 = 1;
du = 0.2;

c1 = sqrt(gamma*p1/rho1);
u1 = M1*c1;

gamma1 = gamma-1;
C  = (M1+1/(gamma*M1))^2/(0.5*M1^2+1/gamma1);

coef = [C/2-1,C/gamma1-2/gamma,-1/gamma^2];
M2 = sqrt(roots(coef));

[bid,indx] = min(abs(M2-M1));
abs(M2(indx)-M1)<1e-10 || error("inconsistent result");

M2(indx) = [];

u2 = u1*sqrt((0.5+1/M1^2/gamma1)/(0.5+1/M2^2/gamma1));
c2 = u2/M2;
rho2 = u1*rho1/u2;
p2 = c2^2*rho2/gamma;

U1 = [rho1,u1+du,0,p1];
U2 = [rho2,u2+du,0,p2];

## dfx = exp(-((x-Lx/2)/sigma).^2);
dfx = (1+sign(x-Lx/5))/2;
Uini(1:nnod,:) = U1(ones(nnod,1),:);
dw = U2-U1;
Uini(1:nnod,:) = Uini(1:nnod,:) + dfx*dw;

asave("gfmovshock.ini.tmp",Uini);
