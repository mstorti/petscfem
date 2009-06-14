## __INSERT_LICENSE__
## $Id: mkstrip2d.m,v 1.1 2003/01/09 02:37:45 mstorti Exp $

source("data.m.tmp");

fricref = aload("friction.step-1.ref");
fric = aload("friction.step-1.tmp");

tol = 1e-10;
erro = merr(fricref-fric);
printf("Test OK ? %d (error %g, tol %g)\n",
       erro<tol,erro,tol);
