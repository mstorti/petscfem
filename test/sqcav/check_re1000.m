##__INSERT_LICENSE__
## $Id: check_re1000.m,v 1.4 2003/01/08 15:49:05 mstorti Exp $
source("data.m.tmp");
ghia;

load -force sqcav.ny.tmp
u=aload("sqcav.weak_form_0.tmp");
uwf1=aload("sqcav.weak_form_1.tmp");
uref = aload("u_re1000.dat");

tol=1e-8;
erro = merr(u(ny,:)-uref);
printf("Weak form 0. error < tol OK ? %d (error = %g, tol = %g)\n", \
       erro<tol,erro,tol);
erro = merr(uwf1(ny,:)-uref);
printf("Weak form 1. error < tol OK ? %d (error = %g, tol = %g)\n", \
       erro<tol,erro,tol);
