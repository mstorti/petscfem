##__INSERT_LICENSE__
## $Id: mksqcav5.m,v 1.2 2003/01/08 15:49:03 mstorti Exp $
N=20; # The number of elements per side
w=zhomo([0 1 0 1],N+1,N+1,[1 5 1 1 5 1]);
[xnod,icone] = pfcm2fem(w);
icone = icone(:,[1 4 3 2]);

asave("sqcav5.nod.tmp",xnod);
asave("sqcav5.con.tmp",icone);

x=xnod(:,1);
y=xnod(:,2);

tol=1e-5;
lid=find(abs(y-1)<tol)';
nlid=length(lid);
b=create_set([find(abs(x)<tol);
              find(abs(x-1)<tol);
              find(abs(y)<tol)]);
b=complement(lid,b);
nb=length(b);

fixa=[lid' ones(nlid,2)*diag([1 1]);
      lid' ones(nlid,2)*diag([2 0]);
      b' ones(nb,2)*diag([1 0]);
      b' ones(nb,2)*diag([2 0])];

asave("sqcav5.fixa.tmp",fixa);
