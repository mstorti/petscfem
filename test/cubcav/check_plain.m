#$Id: check_plain.m,v 1.4 2004/11/11 18:23:08 mstorti Exp $

u0 = aload("cubcav.state.plain_bupl0.tmp");
u1 = aload("cubcav.state.plain_bupl1.tmp");
u2 = aload("cubcav.state.plain_bupl2.tmp");

uref = aload("cubcav.sol");
indx = uref(:,1);
uref = uref(:,2:5);

erro = merr(u1-u0);
tol = 1e-10;
printf("Cubcav/Block uploading error OK? > %d (error %g, tol %g)\n",
       erro<tol,erro,tol);

erro2 = merr(u2-u0);
printf("Cubcav/Block uploading 2 error OK? > %d (error %g, tol %g)\n",
       erro2<tol,erro2,tol);

erro3 = merr(u0(indx,:)-uref);
printf("Cubcav/error wrt.reference solution OK? > %d (error %g, tol %g)\n",
       erro3<tol,erro3,tol);

