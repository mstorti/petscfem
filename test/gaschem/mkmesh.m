##__INSERT_LICENSE__
## $Id: mkmesh.m,v 1.4 2003/11/11 15:41:13 mstorti Exp $
source("data.m.tmp");

## rem(N,2)==0 || warning("N should be even");
Nb_fac=1e-4;

w=zhomo([0 Ly 0 Lx],Ny+1,Nx+1);
[xnod,icone] = pfcm2fem(w);
xnod = xnod(:,[2 1]);

nnod = rows(xnod);

Z1 = ones(nnod,1);
nodfile = [xnod u*Z1 0*Z1 patm+rho_liq*gravity*(Ly-xnod(:,2))];
asave("pool.nod.tmp",nodfile);
asave("pool.con.tmp",icone);

fid = fopen("pool.peri.tmp","w");
fidlat = fopen("pool.left.tmp","w");
for k=1:Ny+1
  node = (k-1)*(Nx+1)+1;
  nodep = k*(Nx+1);
  for l=1:5
    fprintf(fid,"%f %d %d   %f %d %d\n",-1,nodep,l,+1,node,l);
  endfor
  fprintf(fid,"%d %d %f\n",node,1,Nb_fac*Nb);
  fprintf(fid,"%d %d %f\n",node,2,CO);
  fprintf(fid,"%d %d %f\n",node,3,CN);
  fprintf(fid,"%d %d %f\n",node,4,CdO);
  fprintf(fid,"%d %d %f\n",node,5,CdN);
endfor
fclose(fid);
fclose(fidlat);

fid = fopen("pool.fixa.tmp","w");
for k=1:Nx+1
  node = k;
  Nbb = Nb_fac*Nb;
  if abs(xnod(node,1)-x_inject)<L_inject/2
    Nbb = Nb;
  endif
  fprintf(fid,"%d %d %f\n",node,1,Nbb);
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

fidr = fopen("pool.bcc-right.tmp","w");
for k=1:Ny
  node = k*(Nx+1);
  nodep = node + Nx + 1;
  fprintf(fid,"%d %d\n",node,nodep);
endfor
fclose(fidr);

ini = [Nb CO CN CdO CdN];
ini = ini(ones(nnod,1),:);
asave("pool.ini.tmp",ini);
