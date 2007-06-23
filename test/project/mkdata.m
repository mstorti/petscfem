###key mkdata.m
### $Id mstorti-v6-2-8-g27bcb6d Thu Jun 21 13:40:10 2007 -0300$

N = 10;
w = zhomo([0,1,0,1],N+1,N+1);
[xnod,icone] = pfcm2fem(w);
icone = [icone(:,[1,2,3]);
         icone(:,[4,3,1])];

asave("xnod1-2D.tmp",xnod);
asave("icone1-2D.tmp",icone);
u = l2(xnod).^2;
asave("u-2D.tmp",u);

M = 20;
x2 = rand(M,2);
asave("xnod2-2D.tmp",x2);
asave("ua-2D.tmp",l2(x2).^2);
