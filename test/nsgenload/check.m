source("data.m.tmp");

uref = aload(["state." case_name ".ref"]);
u = aload(["state." case_name ".tmp"]);

tol=1e-10;
erro = merr(uref-u);

printf("Error < tol OK ? %d (erro %g, tol %g)\n",erro<tol,erro,tol);
