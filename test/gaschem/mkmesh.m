##__INSERT_LICENSE__
## $Id: mkmesh.m,v 1.6 2003/11/12 14:04:09 mstorti Exp $
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
  fprintf(fidlat,"%d %d %f\n",node,1,Nb_fac*Nb);
  fprintf(fidlat,"%d %d %f\n",node,2,Nb_fac*CO);
  fprintf(fidlat,"%d %d %f\n",node,3,Nb_fac*CN);
  fprintf(fidlat,"%d %d %f\n",node,4,CdO);
  fprintf(fidlat,"%d %d %f\n",node,5,CdN);
endfor
fclose(fid);
fclose(fidlat);

## Imposed values at bottom
## Only impose Nb at inlet
fid = fopen("pool.fixa.tmp","w");
for k=1:Nx+1
  node = k;
  fac = Nb_fac;
  if abs(xnod(node,1)-x_inject)<L_inject/2
    fac = 1;
  endif
  fprintf(fid,"%d %d %f\n",node,1,Nb*fac);
  fprintf(fid,"%d %d %f\n",node,2,CO*fac);
  fprintf(fid,"%d %d %f\n",node,3,CN*fac);
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

ini = [Nb*Nb_fac CO CN CdO CdN];
ini = ini(ones(nnod,1),:);
asave("pool.ini.tmp",ini);
