##__INSERT_LICENSE__
## $Id: check_fs.m,v 1.2 2003/01/08 15:49:04 mstorti Exp $
source("data.m.tmp");

load -force sqcav.ny.tmp
u=aload("sqcav.fractional_step_re400.tmp");
uref=aload("sqcav.fractional_step_re400.ref");
rem(N,2)==0 || error("N should be even!!");
ny=N/2*(N+1)+(1:N+1)';
u=u(ny,:);
erro = merr(u-uref);
tol=1e-10;

printf("Square cavity at Re=400. Error < tol OK ? %d (error = %g, tol = %g)\n", \
       erro<tol,erro,tol);
