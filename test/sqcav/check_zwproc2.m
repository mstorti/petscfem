u = aload("sqcav.zwproc2.tmp");

# check values at the center of the cavity
uref = [-1.76837e-01  1.21197e-03  1.53015e-01];
erro = merr(u(61,:)-uref);
tol=1e-5;

disp("Sq. Cavity with 0weight proc (sbprt=2). ");
printf("error < tol OK ? %d (error = %g, tol = %g)\n", \
       erro<tol,erro,tol);
