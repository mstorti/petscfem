###key circular_mesh.m

parab_param.a = Rext;
parab_param.expo = 4;
parab_param.coe2 = 0.2;

x8 = projectd(Rmean*[-cos(pi/4) sin(pi/4)]',
	      [-1 1]',0.1,"parab");
x3 = projectd(Rmean*[cos(pi/4) sin(pi/4)]',
	      [0 1]',0.1,"parab");
x5 = projectd([L Rext]',
	      [0 1]',0.1,"parab");
refl = diag([1 -1]);

XNOD = [1 Rint*[cos(pi/4)  sin(pi/4)];
	2 Rmean*[cos(pi/4) sin(pi/4)];
	3 x3';
	4 L                Rmean*sin(pi/4);
	5 x5';
	6 Rint*[-cos(pi/4) sin(pi/4)];
	7 Rmean*[-cos(pi/4) sin(pi/4)];
	8 x8';
	9 Rint 0;
	10 Rmean 0;
	11 L 0;
	12 Rint*[cos(pi/4)  -sin(pi/4)];
	13 Rmean*[cos(pi/4) -sin(pi/4)];
	14 (refl*x3)';
	15 L                -Rmean*sin(pi/4);
	16 (refl*x5)';
	17 Rint*[-cos(pi/4) -sin(pi/4)];
	18 Rmean*[-cos(pi/4) -sin(pi/4)];
	19 (refl*x8)';
	20 L  Rext2;
	21 Rmean*cos(pi/4)  Rext2;
	22 Rext2*[-cos(pi/4) sin(pi/4)];
	23 Rext2*[-cos(pi/4) -sin(pi/4)];
	24 Rmean*cos(pi/4)  -Rext2;
	25 L                -Rext2];

mapbou_fun = "parabolic_mapbou";
