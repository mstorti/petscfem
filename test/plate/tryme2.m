Rext = 3;
L = 10;
a = 1.3;

theta = (0:100)'/100*2*pi;
x1=cos(theta);
y1=sin(theta);

x2=Rext*x1;
y2=Rext*y1;

y3=(-100:100)'/100*2*a;
x3 = -a+y3.^2/2/a;

x4=(0:10)'/10*L;
y4 = Rext*ones(size(x4));

plot(x1,y1,x2,y2,x3,y3,x4,y4)
