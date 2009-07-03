## __INSERT_LICENSE__
## $Id: mkstrip2d.m,v 1.1 2003/01/09 02:37:45 mstorti Exp $

source("data.m.tmp");
ok = 1;

## Verify max vel (computation of forces are not involved here)
tol = 1e-5;
u = aload("strip.state.tmp");
u = u(1:N+1,1);
umax = max(u);
umax_ana = rho*L^2*gbody/(8.0*visco);
erro = abs(umax-umax_ana);
printf("umax computed %f, analytic %f, error %f, OK? %d (tol %g)\n",
       umax,umax_ana,erro,erro<tol,tol);
ok &= erro<tol;

## Verify force x
h = L/N;
sigma_ana = rho*gbody*L/2.0;
force = aload("strip3d.force.tmp");
fx = -force(1);
sigma = fx/h^2/2.0;
error = abs(sigma-sigma_ana);
printf("sigma computed %f, analytic %f, error %f, OK? %d (tol %g)\n",
       sigma,sigma_ana,erro,erro<tol,tol);
ok &= erro<tol;

## Moment z
mx = force(6);
mx_ana = sigma_ana*h^2*L;
error = abs(mx-mx_ana);
printf("x-moment computed %f, analytic %f, error %f, OK? %d (tol %g)\n",
       mx,mx_ana,erro,erro<tol,tol);
ok &= erro<tol;

## Null forces and moments
erro = max(abs(force(2:5)));
printf("max null forces computed %f, OK? %d (tol %g)\n",
       erro,erro<tol,tol);
ok &= erro<tol;

printf("Test OK %d\n",ok);
