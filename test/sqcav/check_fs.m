##__INSERT_LICENSE__
## $Id: check_fs.m,v 1.4 2005/02/18 09:18:21 mstorti Exp $
source("data.m.tmp");

load -force sqcav.ny.tmp
u=aload("sqcav.fractional_step_re1000.tmp");
uref=aload("sqcav.fractional_step_re1000.ref");
rem(N,2)==0 || error("N should be even!!");
ny=N/2*(N+1)+(1:N+1)';
u=u(ny,:);
erro = merr(u-uref);
tol=1e-5;

printf("Square cavity at Re=1000. Error < tol OK ? %d (error = %g, tol = %g)\n", \
       erro<tol,erro,tol);
