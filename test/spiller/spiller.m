###key spiller.m

C=0.053547;			# Constants of spiller expr.
E = 1.85;
H1 = 66.50;			# vertical position of spiller top position
L1 = 26.48;			# Distance from top of spiller to start
				# of flat bottom
L2 = 50;			# Flat bottom length
h1 = 2;				# water height at top of spiller

H2 = H1-C*L1^E;			# Height of spiller w.r.t. flat bottom

## generate a series aof points on the free surface by interpolation of
## a spline parallel to the spiller curve near he top of the spiller
## and almost constant far from the spiller

xfs
xfs = [0 0;
       0.1*L1 0;

XNOD = [1 0 H1;
	2 L1 0;
	3 L1+L2;
	4 H1+h1;
	5 
	2 Rmean*[cos(pi/4) sin(pi/4)];
	3 x3';
	4 L                Rmean*sin(pi/4);
	5 x5';
	6 Rint*[-cos(pi/4) sin(pi/4)];
	7 Rmean*[-cos(pi/4) sin(pi/4)];
