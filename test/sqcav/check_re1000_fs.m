source("data.m.tmp");

u=aload("sqcav.fractional_step.tmp");
load -force u_re1000_fs.ref

tol=1e-8;
erro = merr(u(ny,1)-uu);
printf("Fractional step error < tol OK ? %d (error = %g, tol = %g)\n", \
       erro<tol,erro,tol);
