## Copyright (C) 2003 Mario A. Storti
##
## This file is part of PETSc-FEM.
##__INSERT_LICENSE__
## $Id: spiller.m,v 1.7 2003/03/17 21:54:37 mstorti Exp $

## Author: Mario Storti
## Keywords: spiller, mesh
global spiller_data

source("data.m.tmp");

## H = bottom height
## h = water height
## y = H + h position of free surface


H2 = H1-C*L1^E;			# Height of spiller w.r.t. flat bottom
h2 = y2-H2;			# water height at outlet
s.C = C;
s.E = E;
s.H1 = H1;
s.L1 = L1;
s.L2 = L2;
s.h1 = h1;
s.h2 = h2;

spiller_data = s;
spiller_data.H2 = H2;

## generate a series aof points on the free surface by interpolation of
## a spline parallel to the spiller curve near he top of the spiller
## and almost constant far from the spiller

L=L1+L2;
xfs = [0      H1;
       0.1*L1 H1;
       0.2*L1 H1;
       L-0.7*L y2;
       L-0.6*L y2;
       L-0.4*L y2;
       L-0.2*L y2;
       L        y2];

for j=1:3
  xpro = projectd(xfs(j,:)',[0;1],1,"spiller_eq");
  xpro(2) = xpro(2) + h1;
  xfs(j,:) = xpro';
endfor
spiller_data.xfs = xfs;

x5 = projectd([L1;H2],[1;3],10,"fs_eq");

XNOD = [1 0 H1;
	2 L1 H2;
	3 L1+L2 H2;
	4 0 H1+h1;
	5 x5';
	6 L y2];
XNOD = XNOD(:,2:3);

if 0
  xx=(0:100)'/100*L;
  yy=spline(xfs(:,1),xfs(:,2),xx);

  y=spiller_fun(xx);
  plot(xx,yy,xfs(:,1),xfs(:,2),'o',xx,y,XNOD(:,1),XNOD(:,2),'og');
endif

ICONE = [1 2 5 4;
	 2 3 6 5];

H = [1 4 Ny 1 4 1;
     1 2 round(0.6*Nx) 1 0 4;
     2 3 round(0.4*Nx) 1 0 2];

[xnod,icone,mesh] = mesher(XNOD,ICONE,H,"spiller_mapbou");

asave("spiller.nod.tmp",xnod);
asave("spiller.con.tmp",icone);

fs = mesher_bound(mesh,[6 5 4]);
bottom = mesher_bound(mesh,[1 2 3]);
inlet = mesher_bound(mesh,[4 1]);
outlet = mesher_bound(mesh,[3 6]);

## Bottom u=v=0 
fid = fopen("spiller.fixa_bot.tmp","w");
nbot = length(bottom);
for k=bottom(1:nbot)'
  fprintf(fid,"%d %d %f\n",k,1,0.);
  fprintf(fid,"%d %d %f\n",k,2,0.);
endfor
fclose(fid);

## Inlet u=uin, v=0
fid = fopen("spiller.fixa_in.tmp","w");
for k=inlet(1:length(inlet)-1)'
  fprintf(fid,"%d %d %f\n",k,1,uin);
  fprintf(fid,"%d %d %f\n",k,2,0.);
endfor
fclose(fid);

## Outlet p=0., v=0
fid = fopen("spiller.fixa_out.tmp","w");
for k=outlet(2:length(outlet))'
  fprintf(fid,"%d %d %f\n",k,2,0.);
  fprintf(fid,"%d %d %f\n",k,3,0.);
endfor
fclose(fid);

## Compute normals to FS
nfs = length(fs);
normal = xnod(fs(3:nfs),:) - xnod(fs(1:nfs-2),:);
normal = leftscal(1./l2(normal),normal);
normal = [-normal(:,2) normal(:,1)];

## SF  v.n=0  !!NO
fid = fopen("spiller.slip.tmp","w");
for j=2:length(fs)-1
  k= fs(j);
  fprintf(fid,"%f    %d %d   %f  %d %d\n",normal(j-1,1),k,1,normal(j-1,2),k,2);
endfor
fclose(fid);

## SF  p = p_atm = 0.
patm = 0.;
fid = fopen("spiller.patm.tmp","w");
for j=2:length(fs)-1
  k= fs(j);
  fprintf(fid,"%d %d %f\n",k,3,patm);
endfor
fclose(fid);

nnod = rows(xnod);
uini = [1 0 0];
uini = uini(ones(nnod,1),:);
asave("spiller.ini.tmp",uini);
