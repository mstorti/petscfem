wref=[0    0    0    0    0    1    4
      0    0   -2    0    4    0   -8
      0    0    0    0    0   -1    4
      0   -2    0    4    0    0   -8
      4    0    0   -8   -8    0   16
      0    2    0    4    0    0   -8
      0    0    0    0    0   -1    4
      0    0    2    0    4    0   -8
      0    0    0    0    0    1    4]/4;

w = aload("tryme9.output.tmp");

erro = merr(wref-w);
tol = 1e-10;
printf("2D 3x3 cloud  OK? %d, erro=%g, tol=%g\n",erro<tol,erro,tol);
