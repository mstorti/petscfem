## Copyright (C) 2003 Mario A. Storti
##
## This file is part of PETSc-FEM.
##__INSERT_LICENSE__
## $Id: mkwave2.m,v 1.2 2003/05/11 17:13:51 mstorti Exp $

## Author: Mario Storti
## Keywords: wave, mesh

source("data.m.tmp");
w = zhomo([0 Lx 0 h],Nx+1,Ny+1,[1 0 1 1 yratio 1]);
[xnod,icone] = pfcm2fem(w);
icone = icone(:,[1 4 3 2]);

## This is deactivated, now we add a perturbation to the bottom
if 0 && !initia
  x = xnod(:,1);
  y = xnod(:,2);
  xnod(:,2) = y + eta0 * sin(2*pi/Lx*x) .* (y/h);
endif

## Add parabolic bump to the bottom
x = xnod(:,1);
y = xnod(:,2);
xbini = (Lx-L_bump)/2;		# start of bump
xbend = xbini+L_bump;		# end of bump
ybump = choose(x<xbini,0*y,x>xbend,0*y,4*t*(x-xbini).*(xbend-x)/L_bump^2);
xnod(:,2) = y + ybump .* (1-y/h);

asave("wave.nod.tmp",xnod);
asave("wave.con.tmp",icone);

## Initial state
nnod = rows(xnod);
ini = [uini 0 0];
ini = ini(ones(nnod,1),:);
asave("wave.ini.tmp",ini);

## No-slip at bottom + slip at top 
fid  = fopen("wave.fixa_bot.tmp","w");
fid2 = fopen("wave.fixa_top.tmp","w");
fid3 = fopen("wave.patm.tmp","w");
fid4 = fopen("wave.mmv_top.tmp","w");
fid5  = fopen("wave_mmv.fixa_bot.tmp","w");
for k=1:Nx+1
  node = (k-1)*(Ny+1)+1;
  fprintf(fid,"%d %d %f\n",node,1,1.);
  fprintf(fid,"%d %d %f\n",node,2,1.);
  fprintf(fid5,"%d %d %f\n",node,1,0.);
  fprintf(fid5,"%d %d %f\n",node,2,0.);
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
fclose(fid5);

## Periodic
fid = fopen("wave.peri.tmp","w");
fid2 = fopen("wave.mmv_peri.tmp","w");
fid3 = fopen("wave.in.tmp","w");
fid4 = fopen("wave.out.tmp","w");
for k=1:Ny+1
  node_right = k;
  node_left = node_right + Nx*(Ny+1);
  for j=1:3
    fprintf(fid,"%f %d %d   %f %d %d\n",-1.0,node_left,j,1.0,node_right,j);
  endfor
  for j=1:2
    fprintf(fid2,"%f %d %d   %f %d %d\n",-1.0,node_left,j,1.0,node_right,j);
  endfor

  fprintf(fid3,"%d %d   %f\n",node_left,1,uini);
  fprintf(fid3,"%d %d   %f\n",node_left,2,0.);

  fprintf(fid4,"%d %d   %f\n",node_right,3,1.);

endfor
fclose(fid);
fclose(fid2);
fclose(fid3);
fclose(fid4);

fs = (Nx+1:-1:1)'*(Ny+1);
asave("wave.nod_fs.tmp",fs);
nfs = Nx;
normal = [0 1];
normal = normal(ones(Nx+1,1),:);
asave("wave.spines.tmp",normal);
