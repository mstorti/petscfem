u = aload("data.txt.tmp");
u = vec(u');
uu = aload("data.out.tmp");

erro = merr(u-uu);
tol=1e-10;
printf("Read and write vector OK ? %d (error %g, tol %g)\n",
       erro<tol,erro,tol);
