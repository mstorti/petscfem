global Rint Rext L Rmean
Rint=1;
Rext=2;
L = 5;

Ntheta = 40;
Nr=10;
Nx=20;

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
	12 Rint*[cos(pi/4)  -sin(pi/4)];
	13 Rmean*[cos(pi/4) -sin(pi/4)];
	14 Rmean*cos(pi/4)  -Rext;
	15 L                -Rmean*sin(pi/4);
	16 L                -Rext;
	17 Rint*[-cos(pi/4) -sin(pi/4)];
	18 Rmean*[-cos(pi/4) -sin(pi/4)];
	19 Rext*[-cos(pi/4) -sin(pi/4)]	];

XNOD = XNOD(:,2:3);

ICONE = [1 2 7 6;
	 2 3 8 7;
	 2 4 5 3;
	 9 10 2 1;
	 10 11 4 2;
	 12 17 18 13;
	 13 18 19 14;
	 12 13 10 9;
	 13 15 11 10;
	 14 16 15 13;
	 6 7 18 17;
	 7 8 19 18;
	 ];

H = [1 2 Nr/2;
     2 3 Nr/2;
     6 17 Ntheta/4;
     1 6 Ntheta/4;
     17 12 Ntheta/4;
     9 1 Ntheta/8;
     9 12 Ntheta/8;
     10 11 Nx];

[xnod,icone,mesh] = mesher(XNOD,ICONE,H);
external = mesher_bound(mesh,[5 3 8 19 14 16]);
outlet = mesher_bound(mesh,[16 15 11 4 5]);
skin = mesher_bound(mesh,[9 1 6 17 12 9]);
