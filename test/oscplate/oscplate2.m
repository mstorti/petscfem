##__INSERT_LICENSE__
## $Id: oscplate2.m,v 1.1 2003/02/17 12:39:28 mstorti Exp $
source("~/.octaverc");
U=aload("oscsome2.sal");;

printf("size is OK? %d\n",length(U)==600);
maxu=max(abs(U(:,2)));
maxp=max(abs(U(:,4)));
printf("max(u)<tol OK? %d,   (%g)\n",maxu<1e-8,maxu);
printf("max(p)<tol OK? %d,   (%g)\n",maxp<1e-3,maxp);

printf("\n\n\nv en la pared, primer periodo\n");
U(1:2:32,3)

printf("\n\n\nv en x=0.5, primer periodo\n");
U(2:2:32,3)

printf("\n\n\nv en la pared\n");
U(101:20:250,3)

printf("\n\n\nv en x=0.5\n");
U(102:50:600,3)
