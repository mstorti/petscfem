###key check_test_abso.m

uref = aload("canal1d.h-step10");

u=aload("canal1d.rslt.tmp");
u=reshape(u(:,2),103,50);

erro = merr(u(:,10:10:30)-uref);
tol=1e-10;

printf("Canal1d/sw_abso. Comparison wrt. \
reference sol OK? > %d (erro %g, tol %g)\n",erro<tol,erro,tol);
