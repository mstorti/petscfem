u_lu=aload("save.state.lu.tmp");
u_slu=aload("save.state.direct_superlu.tmp");

tol = 1e-8;

err = merr(u_slu-u_lu);
printf("Direct/SuperLU  OK ? > %d, [error: %g]\n",err<tol,err);
