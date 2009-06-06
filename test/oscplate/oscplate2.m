##__INSERT_LICENSE__
## $Id: oscplate2.m,v 1.1 2003/02/17 12:39:28 mstorti Exp $
source("~/.octaverc");
U=aload("oscsome2.sal");;
load oscplate2-ref.octave

printf("size is OK? %d\n",length(U)==600);
maxu=max(abs(U(:,2)));
maxp=max(abs(U(:,4)));
printf("max(u)<tol OK? %d,   (%g)\n",maxu<1e-8,maxu);
printf("max(p)<tol OK? %d,   (%g)\n",maxp<1e-3,maxp);
ok1 = length(U)==600 && maxu<1e-8 && maxp<1e-3;

tol = 1e-5;
ok2 = merr(U(1:2:32,3)-uwall_p1)<tol;
printf("v at wall, first period OK? %d\n",ok2);

ok3 = merr(U(2:2:32,3)-ux05_p1)<tol;
printf("v at x=0.5, first period OK? %d\n",ok3);

ok4 = merr(U(101:20:250,3)-uwall)<tol;
printf("v at wall OK? %d\n",ok4);

ok5 = merr(U(102:50:600,3)-ux05)<tol;
printf("v at x=0.5 OK? %d\n",ok5);

printf("test OK %d\n",ok1 && ok2 && ok3 && ok4 && ok5);

if 0
  uwall_p1 = U(1:2:32,3);
  ux05_p1 = U(2:2:32,3);
  uwall = U(101:20:250,3);
  ux05 = U(102:50:600,3);
  save oscplate2-ref.octave uwall_p1 ux05_p1 uwall ux05
endif
