##__INSERT_LICENSE__
## $Id: sfilter.m,v 1.6 2003/01/10 12:38:54 mstorti Exp $
omega=0.9;
xi=1;
N=128;

w2=0.25;
a2=[1 -(1-w2) 0]';
b2=[w2 0 0]';

x=ones(32,1);
x(1:5)=0;

y = filter(b,a,x);
y2 = filter(b2,a2,x);

t=[1 -1 1];
printf("filter 1 at pi*h: %f\n",(t*b)/(t*a));
printf("filter 2 at pi*h: %f\n",(t*b2)/(t*a2));

f1 = filana(a,b);
f2 = filana(a2,b2);

semilogy(abs([f1(1:N/2) f2(1:N/2)]))
pause
loglog((2:N/2)',abs([f1(2:N/2) f2(2:N/2)]-1))
