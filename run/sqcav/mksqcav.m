N=60; # The number of elements per side
w=zhomo([0 1 0 1],N+1,N+1,[1 5 1 1 5 1]);
[xnod,icone] = pfcm2fem(w);
icone = icone(:,[1 4 3 2]);

asave("sqcavr.nod.tmp",xnod);
asave("sqcavr.con.tmp",icone);

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

asave("sqcavr.fixa.tmp",fixa);
