##__INSERT_LICENSE__
## $Id: mkmesh.m,v 1.1 2003/11/10 20:35:03 mstorti Exp $
source("data.m.tmp");

## rem(N,2)==0 || warning("N should be even");

w=zhomo([0 Ly 0 Lx],Ny+1,Nx+1);
[xnod,icone] = pfcm2fem(w);

asave("pool.nod.tmp",xnod);
asave("pool.con.tmp",icone);

nnod = rows(xnod);

fid = fopen("pool.peri.tmp","w");
for k=1:Ny+1
  node = (k-1)*(Nx+1)+1;
  nodep = k*(Nx+1);
  for l=1:5
    fprintf(fid,"%f %d %d   %f %d %d\n",-1,nodep,l,+1,node,l);
  endfor
endfor
fclose(fid);

