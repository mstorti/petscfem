N = 40;

w = zhomo([0 1 0 1],N+1,N+1,[1 5 1 1 5 1]);
[xnod,icone] = pfcm2fem(w);
icone = icone(:,[1 4 3 2]);
icone_dx = icone(:,[1 2 4 3])-1;

asave("sqcav.nod.tmp",xnod);
asave("sqcav.con.tmp",icone);
asave("sqcav.con0.tmp",icone_dx);

tol = 1e-5;
wall = find(abs(xnod(:,2))<tol | abs(xnod(:,1))<tol | \
	    abs(xnod(:,1)-1)<tol);

pffixa("sqcav.fixa-wall.tmp",wall,1:2);

tol = 1e-5;
top = find(abs(xnod(:,2)-1)<tol);
top = complement(wall,top)';
pffixa("sqcav.fixa-top.tmp",top,1:2,[1 0]);

