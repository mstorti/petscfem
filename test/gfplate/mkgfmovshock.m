## $Id: mkgfmovshock.m,v 1.3 2005/02/08 01:49:15 mstorti Exp $
source("data.m.tmp");

## Nbr of elements along `y' 
Ny = pery;

## Mesh step in `x' direction
hx = Lx/Nx;

## Basic mesh
w = zhomo([0,pery*Lx/Nx,0,Lx],pery+1,Nx+1);
[xnod,icone] = pfcm2fem(w);
xnod = xnod(:,[2 1]);
nnod = size(xnod,1);
x = xnod(1:nnod,1);

## Periodic b.c. along y
m = pery+1;			# nbr of nodes along `y'
pfperi("gfmovshock.peri-y.tmp", \
       pery*(Nx+1)+(1:Nx+1), \
       modulo((0:Nx)+perx,Nx+1)+1,1:4);

## Periodic b.c. along x
pfperi("gfmovshock.peri-x.tmp", \
       (Nx+1)*(1:pery),1+(Nx+1)*(0:pery-1),1:4);

## Save mesh
asave("gfmovshock.nod.tmp",xnod);
asave("gfmovshock.con.tmp",icone);

## Normal to the shock wave
nor = [pery,-perx]';
nor = nor/l2(nor');

## Compute steady shock wave 
M1 = Machin;
p1 = 1;
rho1 = 1;
du = 0.2;

## Values upstream
c1 = sqrt(gamma*p1/rho1);
u1 = M1*c1;

gamma1 = gamma-1;
C  = (M1+1/(gamma*M1))^2/(0.5*M1^2+1/gamma1);

coef = [C/2-1,C/gamma1-2/gamma,-1/gamma^2];
M2 = sqrt(roots(coef));

[bid,indx] = min(abs(M2-M1));
abs(M2(indx)-M1)<1e-10 || error("inconsistent result");

## M2 is Mach downstream
M2(indx) = [];

## Values downstream
u2 = u1*sqrt((0.5+1/M1^2/gamma1)/(0.5+1/M2^2/gamma1));
c2 = u2/M2;
rho2 = u1*rho1/u2;
p2 = c2^2*rho2/gamma;

## Maximum velocity for computing time step
vmax = max([abs(u1)+c1,abs(u2)+c2]);
Dt = hx/vmax;
fid = fopen("gfmovshock.octave-out.tmp","w");
fprintf(fid,"$Dt = %f;\n",Dt);
fclose(fid);

U1 = [rho1,u1+du*nor',p1];
U2 = [rho2,u2+du*nor',p2];

## Coordinate normal to shock
xx = xnod*nor;

nnod = rows(xnod);
Lnor = Lx*nor(1);		# Domain length along
				# normal direction
dfx = (xx>0.1*Lnor) & (xx<0.4*Lnor);
Uini(1:nnod,:) = U2(ones(nnod,1),:);
dw = U1-U2;
Uini(1:nnod,:) = Uini(1:nnod,:) + dfx*dw;

asave("gfmovshock.ini.tmp",Uini);
