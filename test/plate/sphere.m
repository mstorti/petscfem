##__INSERT_LICENSE__
## $Id: sphere.m,v 1.1 2003/01/13 03:30:50 mstorti Exp $
global Rint Rext Rext2 L Rmean

source("data.m.tmp");
rem(Ntheta,8)==0 || error("Ntheta must multiple of 8");
Rmean = (Rint+Rext)/2;

XNOD = [1 Rint*[cos(pi/4)  sin(pi/4)];
	2 Rmean*[cos(pi/4) sin(pi/4)];
	3 Rmean*cos(pi/4)  Rext;
	4 L                Rmean*sin(pi/4);
	5 L                Rext;
	6 Rint*[-cos(pi/4) sin(pi/4)];
	7 Rmean*[-cos(pi/4) sin(pi/4)];
	8 Rext*[-cos(pi/4) sin(pi/4)];
	9 Rint 0;
	10 Rmean 0;
	11 L 0;
	12 -Rint 0;
	13 -Rmean 0;
	14 -Rext 0];

XNOD = XNOD(:,2:3);

ICONE = [1 2 7 6;
	 2 3 8 7;
	 2 4 5 3];

ICONEA=[8 14 13 7;
	7 13 12 6;
	1 9 10 2;
	2 10 11 4];

H = [1 2 Nr/2;
     6 7 Nr/2;
     2 3 Nr/2;
     13 14 Nr/2;
     6 12 Ntheta/8;
     1 6 Ntheta/4;
     1 9 Ntheta/8;
     2 4 Nx;
     10 11 Nx];

[xnod,icone,mesh] = mesher(XNOD,ICONE,H,"mapbous");

[x3d,ic3d] = extrude(xnod,icone,Nphi,1/Nphi);
rho = x3d(:,2);
x   = x3d(:,1);
phi = 2*pi*x3d(:,3);

x3d = [x rho.*cos(phi) rho.*sin(phi)];

nnod = rows(xnod);
[x3d2,ic3d2,iren] = pasten(x3d,ic3d,[(1:nnod)'; nnod*Nphi+(1:nnod)']);

[xnoda,iconea,mesha] = mesher(XNOD,ICONEA,H,"mapbous");

X = [-1 -1;
     +1 -1;
     +1 +1;
     -1 +1]
     
XNODC = 1/sqrt(2)*[0.5*X; X]
ICONE = [1 2 3 4;
	 2 6 7 3;
	 3 7 8 4;
	 4 8 5 1;
	 1 5 6 2];

for k=1:Nr/2
  row = mesher_row(mesha,1,[14 8],k);
  xrow = xnod(row,:);
  
endfor

asave("sphere.nod.tmp",x3d);
asave("sphere.con.tmp",ic3d);

#mesher_row(mesh,1,[1 2],1)
