##__INSERT_LICENSE__
## $Id: checkdlf.m,v 1.2 2003/01/08 15:49:03 mstorti Exp $
source("data.m.tmp");

u=aload(["stream.rot0.tmp"]);
u=reshape(u,Nx+1,nstep);
u=u(1,:)';

uref = aload(["stream.dl_" casen ".ref"]);
erro = merr(u-uref);
tol = 1e-8;
printf("Test OK ? %d. (erro %g, tol %g)\n",erro<tol,erro,tol);
