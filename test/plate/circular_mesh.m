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
	19 Rext*[-cos(pi/4) -sin(pi/4)]	
	20 L  Rext2;
	21 Rmean*cos(pi/4)  Rext2;
	22 Rext2*[-cos(pi/4) sin(pi/4)];
	23 Rext2*[-cos(pi/4) -sin(pi/4)];
	24 Rmean*cos(pi/4)  -Rext2;
	25 L                -Rext2];


mapbou_fun = "circular_mapbou";
