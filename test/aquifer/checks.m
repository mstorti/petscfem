##__INSERT_LICENSE__
## $Id: checks.m,v 1.3 2003/01/08 15:49:03 mstorti Exp $
source("data.m.tmp");

tol=1e-3;

u = aload("chezy_rect.state.tmp");
ref = aload("chezy_rect.state.ref");
erro = merr(u-ref);
printf("Test OK? > %d, [error=%g, tol=%g]\n",erro<tol,erro,tol);

u2 = aload("mann_circ.state.tmp");
ref2 = aload("mann_circ.state.ref");
erro = merr(u2-ref2);
printf("Test OK? > %d, [error=%g, tol=%g]\n",erro<tol,erro,tol);

