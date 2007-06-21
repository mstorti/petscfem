ua = aload("ua-2D.tmp");
ui = aload("ui-2D.tmp");

tol = 8e-3;
erro = merr(ua-ui);
printf("erro %g, tol %g, erro<tol %d\n",erro,tol,erro<tol)
