##__INSERT_LICENSE__
## $Id: oscplate.m,v 1.2 2003/01/08 15:49:03 mstorti Exp $
omega=1000*2*pi;
L=1;
N=10;
nu=100;

T=2*pi/omega;
Dt=T/16.;
t=52*Dt;
x=(0:N)'/N*L;
la1=(1+i)*sqrt(omega/2/nu);
la2=-la1;

phase=exp(i*omega*t);
num=exp(la1*(x-L))-exp(la2*(x-L));
den=exp(-la1*L)-exp(-la2*L);

uex=imag(phase/den*num);

plot(x,uex);
