u1=aload("save.state.iisd_peri.np1.tmp");
u2=aload("save.state.iisd_peri.np2.tmp");

tol = 1e-10;

err = merr(u1-u2);
printf(["IISD on 2 processors with periodic b.c.'s OK ?" \
        " > %d, [error: %g, tol %g]\n"], \
        err<tol,err,tol);
