##__INSERT_LICENSE__
## $Id: sfilter.m,v 1.4.2.1 2003/01/11 13:56:26 mstorti Exp $
omega=0.25;
xi=1;

a = [1/omega^2+xi/omega+1/2;
     -2/omega^2;
     1/omega^2-xi/omega+1/2];

b = [0.5+xi/omega 0 0.5-xi/omega]';

a4 = [a;0];
b4 = [2*b;0]-[0;b];

w2=0.25;
a2=[1 -(1-w2) 0]';
b2=[w2 0 0]';

a3=[1 0 0];
b3=[0 2 -1];

N=128;
f1 = filana(a,b,N);
f2 = filana(a2,b2,N);
f4 = filana(a4,b4,N);

x=ones(N,1);
x(1:3)=0;

yy = filter(b,a,x);
y = filter(b3,a3,yy);

y2 = filter(b2,a2,x);
y3 = filter(b3,a3,y2);

y4 = filter(b4,a4,x);

semilogy(abs([f1 f2 f4]));
pause
indx = (2:N/2)';
loglog(indx,abs([f1(indx) f2(indx) f4(indx)]-1))
pause
plot([y y2])
