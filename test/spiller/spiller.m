###key spiller.m
global spiller_data

## H = bottom height
## h = water height
## y = H + h position of free surface

C=0.053547;			# Constants of spiller expr.
E = 1.85;
H1 = 66.50;			# vertical position of spiller top position
L1 = 26.48;			# Distance from top of spiller to start
				# of flat bottom
L2 = 50;			# Flat bottom length
h1 = 4;				# water height at top of spiller
y2 = 60;			# restitution height

h2 = y2-H2;
s.C = C;
s.E = E;
s.H1 = H1;
s.L1 = L1;
s.L2 = L2;
s.h1 = h1;
s.h2 = h2;

spiller_data = s;
H2 = H1-C*L1^E;			# Height of spiller w.r.t. flat bottom
spiller_data.H2 = H2;

## generate a series aof points on the free surface by interpolation of
## a spline parallel to the spiller curve near he top of the spiller
## and almost constant far from the spiller

L=L1+L2;
xfs = [0      H1;
       0.1*L1 H1;
       0.2*L1 H1;
       L-0.7*L y2;
       L-0.6*L y2;
       L-0.4*L y2;
       L-0.2*L y2;
       L        y2];

for j=1:3
  xpro = projectd(xfs(j,:)',[0;1],1,"spiller_eq");
  xpro(2) = xpro(2) + h1;
  xfs(j,:) = xpro';
endfor
spiller_data.xfs = xfs;

x5 = projectd([L1;H2],[1;1],10,"fs_eq");

XNOD = [1 0 H1;
	2 L1 H2;
	3 L1+L2 H2;
	4 0 H1+h1;
	5 x5';
	6 L y2];

if 1
  xx=(0:100)'/100*L;
  yy=spline(xfs(:,1),xfs(:,2),xx);

  y=spiller_fun(xx);
  plot(xx,yy,xfs(:,1),xfs(:,2),'o',xx,y,XNOD(:,2),XNOD(:,3),'og');
endif
