###key verif.m
### $Id: $

uref = aload("bubble.state.ref");
u = aload("bubble.state.tmp");

tol = 1e-8;
erro = merr(u-uref);
printf("test OK ? %d (erro %f, tol %f)\n",erro<tol,erro,tol);
