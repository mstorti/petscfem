# $Id: check_srfgath.m,v 1.1 2005/02/21 00:35:15 mstorti Exp $

refval=[1 0.12 0.5 0.5 0.12 1 0.35 0.5 0.5 0.35];
vals = aload("cubcav-srfgath.gather.tmp");

erro = merr(refval-vals);
tol = 1e-10;
printf("Srfgath values OK? > %d (error %g, tol %g)\n", \
       erro<tol,erro,tol);
