##__INSERT_LICENSE__
## $Id: check_var_rain.m,v 1.3 2003/01/08 15:49:03 mstorti Exp $
###key check_var_rain.m

source("data.m.tmp");

u=aload("var_rain.some_out.tmp");
u=u(:,2);
nstep=rows(u);
t=(1:nstep)'*Dt;

coef = -rain0/S/eta0;
uex = coef*((t<T0).*t + (t>=T0).*T0.*ones(size(t)));

tol=1e-8;
erro = merr(u-uex);

printf("Variable rain test OK ? %d. (erro %g, tol %g)\n",erro<tol,erro,tol);
