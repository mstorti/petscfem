source("data.m.tmp");
ghia;

load -force sqcav.ny.tmp
u=aload("u.Re1000.tmp");
uref = aload("u_re1000.dat");

tol=1e-8;
erro = merr(u(ny,:)-uref);
printf("error < tol OK ? %d (error = %g, tol = %g)\n",erro<tol,erro,tol);
