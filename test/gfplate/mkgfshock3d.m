## $Id: mkgfshock3d.m,v 1.3 2005/02/27 18:41:51 mstorti Exp $
source("data.m.tmp");

## Tuyere data
plt_file = "~/PETSC/COMP-CORNER/Z23_contour.plt";
z = aload(plt_file);
xtuy = z(:,1);
Rtuy = z(:,2);
x0 = min(z(:,1));
x1 = max(z(:,1));
Ltuy = x1-x0;
abs(Lx-Ltuy)<1e-3 || error("Not matching Ltuy...");

## Make 2D mesh
w = zhomo([0 1 0 Lx ],Nr+1,Nx+1);
[xnod,icone] = pfcm2fem(w);
xnod = xnod(:,[2 1]);
nnod = size(xnod,1);
x = xnod(1:nnod,1);
x = x0+(x/Lx)*Ltuy;

## Scales `y' as a fuction of `x'
xw = x(1:Nx+1,1);
rw = spline(xtuy,Rtuy,xw);
Nx1 = Nx+1;
for k=1:Nx+1
  indx = k+(0:Nr)'*Nx1;
  xnod(indx,2) = xnod(indx,2)*rw(k);
endfor

## Compute derivatives
## in order to impose slip condition
dx = 1e-4*Ltuy;
yyp = spline(xtuy,Rtuy,xw+dx);
yym = spline(xtuy,Rtuy,xw-dx);
dydx = (yyp-yym)/(2*dx);

## Extrude to 3D
[x3d,ico3d] = extrude(xnod,icone,1);
z = x3d(:,1);
rho = x3d(:,2);
theta = dtheta*x3d(:,3);
x3d = [rho.*cos(theta),rho.*sin(theta),z];
nnod3d = rows(x3d);

asave("gfshock3d.nod.tmp",x3d);
asave("gfshock3d.con.tmp",ico3d);
asave("gfshock3d.dx-con.tmp",ico3d(:,[5 6 8 7 1 2 4 3])-1);

## Fixa `y' component of velocity on plane `y=0'
pffixa("gfshock3d.fixa-v.tmp",(1:nnod),3);

## Periodic on `theta=D theta'
pfperi("gfshock3d.peri.tmp",nnod+(1:nnod)',(1:nnod)',[1 4 5]);
n2 = nnod+(1:nnod)';
n1 = (1:nnod)';
cosdt = cos(dtheta);
sindt = sin(dtheta);
pfconstr("gfshock3d.peri-u.tmp", \
	 [n2,n1,n1],[2,2,3], \
	 [-1,cosdt,-sindt]);
pfconstr("gfshock3d.peri-v.tmp", \
	 [n2,n1,n1],[3,2,3], \
	 [-1,sindt,cosdt]);

## rho,u,v at inlet
inlet = 1+(0:Nr)'*Nx1;
pffixa("gfshock3d.fixa-rho-in.tmp", inlet,1,1.0);

## p at inlet (will be scaled by a ramp, though)
inlet = 1+(0:Nr)'*Nx1;
pffixa("gfshock3d.fixa-p-in.tmp", inlet,5,1);

## p at outlet
outlet = Nx1+(0:Nr)'*Nx1;
pffixa("gfshock3d.fixa-p-outlet.tmp", \
       outlet,5,pout0)

## Slip on the axis
pffixa("gfshock3d.axis-slip.tmp", \
       (1:Nx+1)',[2 3]);

## Slip on the wall
wall = Nx1*Nr+(1:Nx1)';
pfconstr("gfshock3d.wall-slip.tmp", \
	 [wall,wall],[4 2], \
	 [-dydx,ones(Nx1,1)]);

Uini = [rhoout0,0,0,0,pout0];
Uini = Uini(ones(nnod3d,1),:);

asave("gfshock3d.ini.tmp",Uini);
asave("gfshock3d.some-nodes.tmp",(1:Nx1)');
