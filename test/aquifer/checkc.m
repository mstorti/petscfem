##__INSERT_LICENSE__
## $Id: checkc.m,v 1.2 2003/01/08 15:49:03 mstorti Exp $
up=aload(["aqui.plain.tmp"]);
uv=aload(["aqui.var_eta0.tmp"]);

erro = merr(up-uv);
tol = 1e-10;

printf("OK ? %d (erro %g, tol %g)\n",erro<tol,erro,tol);
