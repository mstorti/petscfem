## Copyright (C) 2003 Mario A. Storti
##
## This file is part of PETSc-FEM.
##__INSERT_LICENSE__
## $Id: mkwave.m,v 1.1 2003/03/30 15:55:08 mstorti Exp $

## Author: Mario Storti
## Keywords: wave, mesh
source("data.m.tmp");
w = zhomo([0 Lx 0 h],Nx+1,Ny+1,[1 0 1 1 yratio 1]);
[xnod,icone] = pfcm2fem(w);
icone = icone(:,[1 4 3 2]);

asave("wave.nod.tmp",xnod);
asave("wave.con.tmp",icone);

## No-slip at bottom + slip at top
fid = fopen("wave.fixa_bot.tmp","w");
fid2 = fopen("wave.fixa_top.tmp","w");
for k=1:Nx
  node = (k-1)*(Ny+1)+1;
  fprintf(fid,"%d %d %f\n",node,1,0.);
  fprintf(fid,"%d %d %f\n",node,2,0.);
  node = k*(Ny+1);
  fprintf(fid2,"%d %d %f\n",node,2,0.);
endfor
fclose(fid);
fclose(fid2);

## Periodic
fid = fopen("wave.peri.tmp","w");
for k=1:Ny+1
  node_right = k;
  node_left = node_right + Nx*(Ny+1);
  for j=1:3
    fprintf(fid,"%f %d %d   %f %d %d\n",-1.0,node_left,j,1.0,node_right,j);
  endfor
endfor
fclose(fid);
