global Rint Rext L Rmean
Rint=1;
Rext=2;
L = 5;
Ntheta = 20;
Nr=20;
Nx=20;

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
	11 L 0];

XNOD = XNOD(:,2:3);

ICONE = [1 2 7 6;
	 2 3 8 7;
	 2 4 5 3;
	 9 10 2 1;
	 10 11 4 2];

xnod = []; icone=[];

for k=1:rows(ICONE)
  [xx,ii] = isomap(XNOD,ICONE(k,:),20,20);
  nnod = rows(xnod);
  xnod = [xnod; xx];
  icone = [icone; ii+nnod];
endfor
