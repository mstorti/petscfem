##$Id#
###key verif_fd_jac_1.m

uref = aload("save.state.fd_jac_1_ref.tmp");
u = aload("save.state.fd_jac_1.tmp");

erro=merr(u-uref);
tol = 1e-10;
printf("Solution with fd jac. computed OK ? %d\n",erro<tol);
printf("error %g, tol %g\n",erro,tol);
