##__INSERT_LICENSE__
## $Id: check_zwproc.m,v 1.3 2003/01/08 15:49:05 mstorti Exp $
u = aload("sqcav.zwproc.tmp");
uref = aload("sqcav.zwproc_ref.tmp");

# check values at the center of the cavity
erro = merr(u-uref);
tol=1e-10;

printf("Sq. Cavity with 0weight proc. error < tol OK ? %d (error = %g, tol = %g)\n", \
       erro<tol,erro,tol);
