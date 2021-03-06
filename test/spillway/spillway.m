## Copyright (C) 2003 Mario A. Storti
##
## This file is part of PETSc-FEM.
##__INSERT_LICENSE__
## $Id: spillway.m,v 1.10 2003/03/30 15:55:08 mstorti Exp $

## Author: Mario Storti
## Keywords: spillway, mesh
global spillway_data

source("data.m.tmp");
initia = 1;
if !exist("initia");
  initia = str2num(getenv("initia"));
endif
## printf("# initia: %d\n",initia);

## H = bottom height
## h = water height
## y = H + h position of free surface

H2 = H1-C*L1^E;			# Height of spillway w.r.t. flat bottom
h2 = y2-H2;			# water height at outlet
s.C = C;
s.E = E;
s.H1 = H1;
s.L1 = L1;
s.L2 = L2;
s.h1 = h1;
s.h2 = h2;

spillway_data = s;
spillway_data.H2 = H2;

## generate a series aof points on the free surface by interpolation of
## a spline parallel to the spillway curve near he top of the spillway
## and almost constant far from the spillway

L=L1+L2;

xfs = onedstr([1 0 1],10)*L1;
xfs(length(xfs))=[];
xfs = [xfs;
       L1+onedstr([1 0 5],20)*L2];
npc = rows(xfs);

if initia
  xfs1 = [0      H1;
	  0.025*L1 H1;
	  0.05*L1 H1;
	  0.10*L1 H1;
	  0.2*L1 H1;
	  0.4*L1 H1;
	  L-0.7*L y2;
	  L-0.6*L y2;
	  L-0.4*L y2;
	  L-0.2*L y2;
	  L        y2];
  
  for j=1:6
    xpro = projectd(xfs1(j,:)',[0;1],1,"spillway_eq");
    xpro(2) = xpro(2) + h1;
    xfs1(j,:) = xpro';
  endfor
  xfs = [xfs spline(xfs1(:,1),xfs1(:,2),xfs)];
else
  xnod = aload("spillway.nod.tmp");
  xfs_old = aload("spillway.xfs.tmp");

  if fs_relax==0
    asave("spillway.nod.tmp",xnod);
    asave("spillway.xfs.tmp",xfs_old);
    return;
  endif

  fs = aload("spillway.nod_fs.tmp");
  state = aload("spillway.state.tmp");
  vfs = state(fs,1:2);
  xfs_fem = xnod(fs,:);
  clear xnod
  nfs = length(fs);
  normal_e = xfs_fem(2:nfs,:)-xfs_fem(1:nfs-1,:);
  le = l2(normal_e);
  normal_e = leftscal(1./le,normal_e);
  normal_e = [normal_e(:,2) -normal_e(:,1)];
  normal = leftscal(le(1:nfs-2),normal_e(1:nfs-2,:)) \
      +leftscal(le(2:nfs-1),normal_e(2:nfs-1,:));
  normal = leftscal(1./l2(normal),normal);
  vn = sum((vfs(2:nfs-1,:).*normal)')';
  vn = [0;vn;0];
  vn_pc = spline(xfs_fem(:,1),vn,xfs);
  xfs = xfs_old;
  normal_pc = zeros(npc,2);
  normal_pc(1,:) = [0 1];
  normal_pc(npc,:) = [0 1];
  epsil = 1e-6*L;
  for k=2:npc-1
    xx = xfs(k,1);
    fp = spline(xfs(:,1),xfs(:,2),xx+epsil);
    fm = spline(xfs(:,1),xfs(:,2),xx-epsil);
    dfdx = (fp-fm)/(2*epsil);
    normal_pc(k,:) = [-dfdx 1]/sqrt(1+dfdx^2);
  endfor
  xfs(2:npc-1,:) = xfs(2:npc-1,:) + fs_relax \
      * Dt*leftscal(vn_pc(2:npc-1),normal_pc(2:npc-1,:));
  pfstep = getenv("petscfem_step");
  asave(["spillway.xfs_" pfstep ".tmp"],xfs);
  ##  save spillway.tmp xfs vn_pc 
  printf("spillway.m: convergence on free surface control points: %g\n",
	 merr(xfs-xfs_old));
endif
asave("spillway.xfs.tmp",xfs);

spillway_data.xfs = xfs;

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

  y=spillway_fun(xx);
  plot(xx,yy,xfs(:,1),xfs(:,2),'o', \
       xx,y,XNOD(:,1),XNOD(:,2),'og');
endif

ICONE = [1 2 5 4;
	 2 3 6 5];

H = [1 4 Ny 1 4 1;
     1 2 round(0.6*Nx) 1 0 4;
     2 3 round(0.4*Nx) 1 0 2];

[xnod,icone,mesh] = mesher(XNOD,ICONE,H,"spillway_mapbou");

asave("spillway.nod.tmp",xnod);
if ~initia; return; endif
asave("spillway.con.tmp",icone);

## Elements with FS damper 
xele = pfnd2ele(xnod,icone,xnod(:,1));
index = find(xele>x_damp);
asave("spillway.damp_con.tmp",icone(index,:));

fs = mesher_bound(mesh,[6 5 4]);
asave("spillway.nod_fs.tmp",fs);
bottom = mesher_bound(mesh,[1 2 3]);
inlet = mesher_bound(mesh,[4 1]);
outlet = mesher_bound(mesh,[3 6]);

## Spines
nfs = length(fs);
xfs_fem = xnod(fs,:);
normal_e = xfs_fem(2:nfs,:)-xfs_fem(1:nfs-1,:);
le = l2(normal_e);
normal_e = leftscal(1./le,normal_e);
normal_e = [normal_e(:,2) -normal_e(:,1)];
normal = leftscal(le(1:nfs-2),normal_e(1:nfs-2,:)) \
    +leftscal(le(2:nfs-1),normal_e(2:nfs-1,:));
normal = leftscal(1./l2(normal),normal);
normal = [0 1;
	  normal;
	  0 1];
asave("spillway.spines.tmp",normal);

## Bottom u=v=0 
fid = fopen("spillway.fixa_bot.tmp","w");
nbot = length(bottom);
for k=bottom(1:nbot)'
  fprintf(fid,"%d %d %f\n",k,1,0.);
  fprintf(fid,"%d %d %f\n",k,2,0.);
endfor
fclose(fid);

## Inlet u=uin, v=0
fid = fopen("spillway.fixa_in.tmp","w");
fid2 = fopen("spillway.mmv_fixa_in.tmp","w");
for k=inlet(1:length(inlet)-1)'
  fprintf(fid,"%d %d %f\n",k,1,uin);
  fprintf(fid,"%d %d %f\n",k,2,0.);
  fprintf(fid2,"%d %d %f\n",k,1,0.);
  fprintf(fid2,"%d %d %f\n",k,2,0.);
endfor
fclose(fid);
fclose(fid2);

## Compute normals to FS
nfs = length(fs);
normal = xnod(fs(3:nfs),:) - xnod(fs(1:nfs-2),:);
normal = leftscal(1./l2(normal),normal);
normal = [-normal(:,2) normal(:,1)];

## SF  v.n=0  !!NO
fid = fopen("spillway.slip.tmp","w");
for j=2:length(fs)-1
  k= fs(j);
  fprintf(fid,"%f    %d %d   %f  %d %d\n",
	  normal(j-1,1),k,1,normal(j-1,2),k,2);
endfor
fclose(fid);

## SF  p = p_atm = 0.
p_atm_var = length(steady_state_file)>0;
if p_atm_var
  state = aload(steady_state_file);
endif
fid = fopen("spillway.patm.tmp","w");
fid2 = fopen("spillway.mmv_fixa_fs.tmp","w");
for j=2:length(fs)-1
  k= fs(j);
  p_atm_node = patm;
  if p_atm_var; p_atm_node = state(k,3); endif
  fprintf(fid,"%d %d %f\n",k,3,p_atm_node);

  fprintf(fid2,"%d %d %f\n",k,1,1.);
  fprintf(fid2,"%d %d %f\n",k,2,1.);
endfor
fclose(fid);
fclose(fid2);

## Outlet p = hydrostatic
p_corner_node = 0;		# Pressure at the intersection(FS,outlet)
if p_atm_var
  nnod = rows(state);
  p_corner_node = state(nnod,3); # The last node is the node at the corner
endif
fid = fopen("spillway.fixa_out.tmp","w");
fid2 = fopen("spillway.mmv_fixa_out.tmp","w");
fid3 = fopen("spillway.bcconv_out.tmp","w");
rho = 1;
corner_node = outlet(length(outlet));
for j=2:length(outlet)
  k1=outlet(j-1);
  k=outlet(j);
  fprintf(fid,"%d %d %f\n",k,3,-rho*gravity*(xnod(k,2)-xnod(corner_node,2)));
  ## fprintf(fid,"%d %d %f\n",k,3,-gravity*xnod(k,2));
  fprintf(fid2,"%d %d %f\n",k,1,0.);
  fprintf(fid2,"%d %d %f\n",k,2,0.);
  fprintf(fid3,"%d %d\n",k1,k);
endfor
fprintf(fid,"%d %d %f\n",corner_node,2,0.);
fclose(fid);
fclose(fid2);
fclose(fid3);

nnod = rows(xnod);
uini = [uin 0 0];
uini = uini(ones(nnod,1),:);
asave("spillway.ini.tmp",uini);
