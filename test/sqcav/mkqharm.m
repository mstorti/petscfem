##__INSERT_LICENSE__
## $Id: mkqharm.m,v 1.2 2003/01/08 15:49:05 mstorti Exp $
source("data.m.tmp");

rem(N,2)==0 || warning("N should be even");

w=zhomo([-1 1 0 2],N+1,N+1,[1 hratio 1 1 hratio 1]);
[xnod,icone] = pfcm2fem(w);
icone = icone(:,[1 4 3 2]);

asave("qharm.nod.tmp",xnod);
asave("qharm.con.tmp",icone);

x=xnod(:,1);
y=xnod(:,2);

tol=1e-5;
wall=find(abs(y)<tol)';
nwall=length(wall);
b=create_set([find(abs(x)<tol);
              find(abs(x-1)<tol);
              find(abs(y)<tol)]);
b=complement(wall,b);
nb=length(b);

fixa=[wall' ones(nwall,2)*diag([1 1]);
      wall' ones(nwall,2)*diag([2 0]);
      b' ones(nb,2)*diag([1 0]);
      b' ones(nb,2)*diag([2 0])];

fid = fopen("qharm.fixa.tmp","w");
fid>0 || error("couldn't open qharm.fixa.tmp");
for k=wall
  xx = x(k);
  fprintf(fid,"%d %d %f\n",k,1,tanh(5.*xx));
endfor
fclose(fid);
