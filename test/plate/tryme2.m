Rext = 3;
L = 10;
a = 2;
expo = 4;

theta = (0:100)'/100*2*pi;
x1=cos(theta);
y1=sin(theta);

x2=Rext*x1;
y2=Rext*y1;

## xmax = (((1+L/a)*2)^(1/expo))*a;
xmax = 4.2;
y3=(-100:100)'/100*xmax;
x3 = a*(-1+(y3/a).^2/2+0.2*abs(y3/a).^expo);

x4=(0:10)'/10*L;
y4 = Rext*ones(size(x4));

plot(x1,y1,x2,y2,x3,y3,x4,y4)
