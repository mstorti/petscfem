##__INSERT_LICENSE__
## $Id: check.g_body.m,v 1.2 2003/01/08 15:49:04 mstorti Exp $
source("data.m.tmp");

x = aload("sqcav.nod.tmp");
u = aload("sqcav.g_body.tmp");

## The solution should be only the hydrostatic pressure
tol=1e-10;
errou = merr(u(:,1));
errov = merr(u(:,2));
errop = merr(u(:,3)+x(:,2));

printf("G_body test, max err in 'u' < tol OK ? %d (error = %g, tol = %g)\n", \
       errou<tol,errou,tol);

printf("G_body test, max err in 'v' < tol OK ? %d (error = %g, tol = %g)\n", \
       errov<tol,errov,tol);

printf("G_body test, max err in 'p+z' < tol OK ? %d (error = %g, tol = %g)\n", \
       errop<tol,errop,tol);

