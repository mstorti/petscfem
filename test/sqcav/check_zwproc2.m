u = aload("sqcav.zwproc2.tmp");
uref = aload("sqcav.zwproc2_ref.tmp");

# check values at the center of the cavity
erro = merr(u-uref);
tol=1e-5;

disp("Sq. Cavity with 0weight proc (sbprt=2). ");
printf("error < tol OK ? %d (error = %g, tol = %g)\n", \
       erro<tol,erro,tol);
