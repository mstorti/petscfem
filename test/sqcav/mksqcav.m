##__INSERT_LICENSE__
## $Id: mksqcav.m,v 1.5 2003/02/12 00:36:27 mstorti Exp $
source("data.m.tmp");

## rem(N,2)==0 || warning("N should be even");

w=zhomo([0 1 0 1],N+1,N+1,[1 hratio 1 1 hratio 1]);
[xnod,icone] = pfcm2fem(w);
icone = icone(:,[1 4 3 2]);

if use_triangles
  icone = [icone(:,[1 2 3]);
	   icone(:,[3 4 1])];
endif

asave("sqcav.nod.tmp",xnod);
asave("sqcav.con.tmp",icone);

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

if g_body 
  utop = 0;
else
  utop = 1;
endif

fixa=[lid' ones(nlid,2)*diag([1 utop]);
      lid' ones(nlid,2)*diag([2 0]);
      b' ones(nb,2)*diag([1 0]);
      b' ones(nb,2)*diag([2 0])];

asave("sqcav.fixa.tmp",fixa);

## y coordinates of nodes on the centerline
ny = (N/2)*(N+1)+(1:N+1)';
yh = xnod(ny,2);

save sqcav.ny.tmp  ny yh
