###key check_test_abso.m

uref = aload("state_h_step10.dat");

u=aload("canal1d.rslt.tmp");
u=reshape(u(:,2),201,51);

erro = merr(u(:,10:10:30)-uref);
tol=1e-9;

printf("Canal1d/sw_abso. Comparison wrt. \
reference sol OK? > %d (erro %g, tol %g)\n",erro<tol,erro,tol);
