##__INSERT_LICENSE__
## $Id: check_peri.m,v 1.2 2003/01/08 15:49:04 mstorti Exp $
u1=aload("save.state.iisd_peri.np1.tmp");
u2=aload("save.state.iisd_peri.np2.tmp");

tol = 1e-10;

err = merr(u1-u2);
printf(["IISD on 2 processors with periodic b.c.'s OK ?" \
        " > %d, [error: %g, tol %g]\n"], \
        err<tol,err,tol);
