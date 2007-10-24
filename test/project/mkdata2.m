###key mkdata.m
### $Id merge-with-petsc-233-50-g0ace95e Fri Oct 19 17:49:52 2007 -0300$

N = 20;
w = zhomo([0,1,0,1],N+1,N+1);
[xnod,icone] = pfcm2fem(w);
[x3d,ico3d] = extrude(xnod,icone,N,1/N);
asave("xnod1.tmp",x3d);
asave("icone1.tmp",ico3d);
system("make convert");

u = l2(x3d).^2;
asave("u.tmp",u);

M = 1.0*5*N^3;
x2 = rand(M,3);
asave("xnod2.tmp",x2);

asave("ua.tmp",l2(x2).^2);
