##__INSERT_LICENSE__
## $Id: checku.m,v 1.2 2003/01/08 15:49:03 mstorti Exp $
source("data.m.tmp");

alpha=K/S;

u=aload(["aqui." casen ".tmp"]);
x=aload("aqui.nod.tmp");

ix=(1:2:2*Nx+2)';
x=x(ix,1);
u=u(ix);

uex = erfc((1-x)/sqrt(4*alpha*t));
## plot(x,[u uex])

erro = max(abs(u-uex));
tol = 1e-2;
printf("test OK ? %d [max err %g, tol %g]\n",erro<tol,erro,tol);
