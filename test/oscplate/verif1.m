###key verif_oscplate1.m
### $Id: $

uref = aload("oscplate1.state.ref");
u = aload("oscplate1.state.tmp");

tol = 1e-5;
erro = merr(uref-u);
printf("|u-uref|<tol OK? %d (error %f, tol %f)\nOK %d\n",
       erro<tol,erro,tol,erro<tol);
