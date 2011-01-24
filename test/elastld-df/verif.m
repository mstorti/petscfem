###key verif.m
### $Id: $

uref = aload("elastld.state40.ref");
u = aload("elastld.state.tmp");

tol = 1e-8;
erro = merr(u-uref);
printf("test OK ? %d (erro %f, tol %f)\n",erro<tol,erro,tol);
