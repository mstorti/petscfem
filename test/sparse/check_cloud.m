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
w1=w(1:9,:);
w2=w(9+(1:9),:);

erro = merr(wref-w1);
tol = 1e-10;
printf("2D 3x3 cloud  OK? %d, erro=%g, tol=%g\n",erro<tol,erro,tol);

erro = merr(wref-w2);
tol = 1e-10;
printf("2D 3x3 cloud (x0!=0) OK? %d, erro=%g, tol=%g\n",erro<tol,erro,tol);
