u = aload("sqcav.zwproc.tmp");

# check values at the center of the cavity
uref = [-1.7731e-01  5.6593e-04  3.3084e-01 ];
erro = merr(u(61,:)-uref);
tol=1e-5;

printf("Sq. Cavity with 0weight proc. error < tol OK ? %d (error = %g, tol = %g)\n", \
       erro<tol,erro,tol);
