##__INSERT_LICENSE__
## $Id: check_g3d.m,v 1.1 2003/02/24 00:14:24 mstorti Exp $

uref = aload("qharm.state_g3d_ref.tmp");
u0 = aload("qharm.state_g3d_0.tmp");
u45 = aload("qharm.state_g3d_45.tmp");
u90 = aload("qharm.state_g3d_90.tmp");

tol=1e-8;
erro  = merr(u0-uref);
printf("angle 0deg, OK ? %d, (error %g, tol %g)\n",erro<tol,erro,tol);

erro  = merr(u45-uref);
printf("angle 45deg, OK ? %d, (error %g, tol %g)\n",erro<tol,erro,tol);

erro  = merr(u90-uref);
printf("angle 90deg, OK ? %d, (error %g, tol %g)\n",erro<tol,erro,tol);
