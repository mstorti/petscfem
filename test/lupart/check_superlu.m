u_lu=aload("save.state.lu.tmp");
u_iisd1=aload("save.state.superlu.np1.tmp");
u_iisd2=aload("save.state.superlu.np2.tmp");
u_iisd2r=aload("save.state.superlu.np2.rand.tmp");
u_iisd_cgs=aload("save.state.superlu.cgs.np2.tmp");

tol = 1e-8;

err = merr(u_iisd1-u_lu);
printf("IISD/SuperLU on 1 processors OK ? > %d, [error: %g]\n",err<tol,err);

err = merr(u_iisd2-u_lu);
printf("IISD/SuperLU on 2 processors OK ? > %d, [error: %g]\n",err<tol,err);

err = merr(u_iisd2r-u_lu);
printf("IISD/SuperLU on 2 processors with rand part. OK ? > %d, [error: %g]\n",err<tol,err);

err = merr(u_iisd_cgs-u_lu);
printf("IISD/SuperLU on 2 processors with CGS OK ? > %d, [error: %g]\n",err<tol,err);
