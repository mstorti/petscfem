#$Id: check_plain.m,v 1.1 2003/12/08 12:33:08 mstorti Exp $

u0 = aload("cubcav.state.plain_bupl0.tmp");
u1 = aload("cubcav.state.plain_bupl1.tmp");

erro = merr(u1-u0);
tol = 1e-10;
printf("Cubcav/Block uploading error OK? > %d (error %g, tol %g)\n",
       erro<tol,erro,tol);
