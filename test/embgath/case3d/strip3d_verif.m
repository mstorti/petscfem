## __INSERT_LICENSE__
## $Id: mkstrip2d.m,v 1.1 2003/01/09 02:37:45 mstorti Exp $

source("data.m.tmp");
ok = 1;

tol = 1e-5;
u = aload("strip.state.tmp");
u = u(1:N+1,1);
umax = max(u);
umax_ana = L^2*gbody/(8.0*visco);
erro = abs(umax-umax_ana);
printf("umax computed %f, analytic %f, error %f, OK? %d (tol %g)\n",
       umax,umax_ana,erro,erro<tol,tol);
ok |= erro<tol;

h = L/N;
sigma_ana = gbody*L;
force = aload("strip3d.force.tmp");
fx = -force(1);
sigma = fx/h^2;
error = abs(sigma-sigma_ana);
printf("sigma computed %f, analytic %f, error %f, OK? %d (tol %g)\n",
       sigma,sigma_ana,erro,erro<tol,tol);
ok |= erro<tol;

printf("Test OK %d\n",ok);
