ua = aload("ua.tmp");
ui = aload("ui.tmp");

tol = 3e-3;
erro = merr(ua-ui);
printf("erro %g, tol %g, erro<tol %d\n",erro,tol,erro<tol)
