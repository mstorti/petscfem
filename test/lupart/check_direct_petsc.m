u_lu=aload("save.state.lu.tmp");
u_petsc=aload("save.state.direct_petsc.tmp");

tol = 1e-7;

err = merr(u_petsc-u_lu);
printf("Direct/PETSc  OK ? > %d, [error: %g]\n",err<tol,err);
