## $Id: mkgfmovshock.m,v 1.6 2005/02/08 03:00:57 mstorti Exp $
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
fprintf(fid,"$norx = %f;\n",nor(1));
fprintf(fid,"$nory = %f;\n",nor(2));
fprintf(fid,"1;\n");
fclose(fid);

U1 = [rho1,u1+du*nor',p1];
U2 = [rho2,u2+du*nor',p2];

nnod = rows(xnod);
Lnor = Lx*nor(1);		# Domain length along
				# normal direction
asave("gfmovshock.dx-nod.tmp",xnod);

nnod2 = nnod;
xx = xnod*nor;

## Aborbing nodes at outlet
tol = 1e-7;
nodes_out = find(xx>Lnor-tol);	# Nodes at outlet
nficout = length(nodes_out);
ficout = nnod+(1:nficout)';
nnod2 = nnod2 + nficout;
xnod = [xnod;
	xnod(nodes_out,:)];
asave("gfmovshock.abso-con-out.tmp", \
      [nodes_out,ficout]);

## Aborbing nodes at outlet
tol = 1e-7;
nodes_out = find(xx>Lnor-tol);	# Nodes at outlet
nficout = length(nodes_out);
ficout = nnod2+(1:nficout)';
nnod2 = nnod2 + nficout;
xnod = [xnod;
	xnod(nodes_out,:)];
asave("gfmovshock.abso-con-out.tmp", \
      [nodes_out,ficout,ficout]);

## Aborbing nodes at inlet
nodes_in = find(xx<+tol);	# Nodes at inlet
nficin = length(nodes_in);
ficin = nnod2 + (1:nficin)';
nnod2 = nnod2 + nficin;
xnod = [xnod;
	xnod(nodes_in,:)];
asave("gfmovshock.abso-con-in.tmp", \
      [nodes_in,ficin,ficin]);

## Coordinate normal to shock
xx = xnod*nor;
dfx = (xx>0.3*Lnor);
Uini = U1(ones(nnod2,1),:);
dw = U2-U1;
Uini = Uini + dfx*dw;

Uini(ficout,:) = 0.;
Uini(ficin,:) = 0.;

asave("gfmovshock.ini.tmp",Uini);

## Save mesh
asave("gfmovshock.nod.tmp",xnod);
asave("gfmovshock.con.tmp",icone);
