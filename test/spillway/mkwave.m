## Copyright (C) 2003 Mario A. Storti
##
## This file is part of PETSc-FEM.
##__INSERT_LICENSE__
## $Id: mkwave.m,v 1.3 2003/03/31 00:01:08 mstorti Exp $

## Author: Mario Storti
## Keywords: wave, mesh
source("data.m.tmp");
w = zhomo([0 Lx 0 h],Nx+1,Ny+1,[1 0 1 1 yratio 1]);
[xnod,icone] = pfcm2fem(w);
icone = icone(:,[1 4 3 2]);

if !initia
  x = xnod(:,1);
  y = xnod(:,2);
  xnod(:,2) = y + eta0 * sin(2*pi/Lx*x) .* (y/h);
endif

asave("wave.nod.tmp",xnod);
asave("wave.con.tmp",icone);

## No-slip at bottom + slip at top 
fid  = fopen("wave.fixa_bot.tmp","w");
fid2 = fopen("wave.fixa_top.tmp","w");
fid3 = fopen("wave.patm.tmp","w");
fid4 = fopen("wave.mmv_top.tmp","w");
for k=1:Nx
  node = (k-1)*(Ny+1)+1;
  fprintf(fid,"%d %d %f\n",node,1,0.);
  fprintf(fid,"%d %d %f\n",node,2,0.);
  node = k*(Ny+1);
  fprintf(fid2,"%d %d %f\n",node,2,0.);
  fprintf(fid3,"%d %d %f\n",node,3,0.);
  fprintf(fid4,"%d %d %f\n",node,1,1.);
  fprintf(fid4,"%d %d %f\n",node,2,1.);
endfor
fclose(fid);
fclose(fid2);
fclose(fid3);
fclose(fid4);

## Periodic
fid = fopen("wave.peri.tmp","w");
fid2 = fopen("wave.mmv_peri.tmp","w");
for k=1:Ny+1
  node_right = k;
  node_left = node_right + Nx*(Ny+1);
  for j=1:3
    fprintf(fid,"%f %d %d   %f %d %d\n",-1.0,node_left,j,1.0,node_right,j);
  endfor
  for j=1:2
    fprintf(fid2,"%f %d %d   %f %d %d\n",-1.0,node_left,j,1.0,node_right,j);
  endfor
endfor
fclose(fid);
fclose(fid2);

fs = (Nx:-1:1)'*(Ny+1);
asave("wave.nod_fs.tmp",fs);
nfs = Nx;
normal = [0 1];
normal = normal(ones(Nx,1),:);
asave("wave.spines.tmp",normal);
