##__INSERT_LICENSE__
## $Id: mkmesh.m,v 1.2 2003/11/10 21:29:50 mstorti Exp $
source("data.m.tmp");

## rem(N,2)==0 || warning("N should be even");

w=zhomo([0 Ly 0 Lx],Ny+1,Nx+1);
[xnod,icone] = pfcm2fem(w);
xnod = xnod(:,[2 1]);

nodfile = [xnod 0*xnod rho_liq*gravity*(Ly-xnod(:,2))];
asave("pool.nod.tmp",nodfile);
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

fid = fopen("pool.fixa.tmp","w");
for k=1:Nx+1
  node = k;
  fprintf(fid,"%d %d %f\n",node,1,Nb);
  fprintf(fid,"%d %d %f\n",node,2,CO);
  fprintf(fid,"%d %d %f\n",node,3,CN);
  fprintf(fid,"%d %d %f\n",node,4,CdO);
  fprintf(fid,"%d %d %f\n",node,5,CdN);
endfor
fclose(fid);

fid = fopen("pool.bcc-up.tmp","w");
for k=1:Nx
  node = Ny*(Nx+1)+k;
  fprintf(fid,"%d %d\n",node+1,node);
endfor
fclose(fid);

ini = [Nb CO CN CdO CdN];
ini = ini(ones(nnod,1),:);
asave("pool.ini.tmp",ini);
