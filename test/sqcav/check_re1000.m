source("data.m.tmp");
ghia;

load -force sqcav.ny.tmp
u=aload("u.Re1000.wf0.tmp");
uwf1=aload("u.Re1000.wf1.tmp");
uref = aload("u_re1000.dat");

tol=1e-8;
erro = merr(u(ny,:)-uref);
printf("Weak form 0. error < tol OK ? %d (error = %g, tol = %g)\n", \
       erro<tol,erro,tol);
erro = merr(uwf1(ny,:)-uref);
printf("Weak form 1. error < tol OK ? %d (error = %g, tol = %g)\n", \
       erro<tol,erro,tol);
