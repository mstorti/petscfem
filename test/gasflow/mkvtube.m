##__INSERT_LICENSE__
## $Id: mkvtube.m,v 1.6 2003/01/20 18:51:13 mstorti Exp $
source("data.m.tmp");

XNOD = [1 0 Rin;
	2 L0 Rin;
	3 L0 R0;
	4 0 R0];

XNOD = XNOD(:,2:3);

ICONE = [1 2 3 4];

H = [1 2 Nz 1 z_ratio 1;
     2 3 Nr r_ratio 0 1];

[xnod,icone,mesh] = mesher(XNOD,ICONE,H,"mapbouvt");
nnod = rows(xnod);

hr = R0/Nr;			# r mesh size
hz = L0/Nz;			# z msh size
h = sqrt(hr*hz);		# typical mesh size
r_av = sqrt(R0*Rin);		# typical radius
Dphi = h/r_av/(2*pi);
if axisymm
  [x3d,ic3d] = extrude(xnod,icone,1,Dphi);
else
  [x3d,ic3d] = extrude(xnod,icone,Nphi,2*pi/Nphi);
endif  
z = x3d(:,1);
rho = x3d(:,2);
phi = x3d(:,3);

x3d = [rho.*cos(phi) rho.*sin(phi) z];

## Paste faces at phi=0, 2*pi
if !axisymm
  iren = (1:rows(x3d))';
  iren(Nphi*nnod+(1:nnod)) = (1:nnod)';
  for k=1:8
    ic3d(:,k) = iren(ic3d(:,k));
  endfor
  x3d(Nphi*nnod+(1:nnod),:) = [];
endif

asave("vtube.nod.tmp",x3d);
asave("vtube.con.tmp",ic3d);

fid = fopen("vtube.fixa.tmp","w");
tol = 1e-7;
#wall = find((abs(rho-R0)<tol || z<tol || abs(z-L0)<tol)
ndim = 3;
u_dof = 1;
if compressible
  u_dof = 2;
endif
p_dof = u_dof+ndim;

done = 0;
nn = rows(x3d);			# total number of nodes
if axisymm; nn=nn/2; endif	# number of nodes in the first layer

if !compressible
  p_h=p_c=0;
endif

n_in=n_h=n_c=n_wall=0;
for k=1:nn
  if k/rows(x3d) > done+0.1;
    done = done+0.1;
    printf("%3d%% done\n",round(100*done));
  endif
  is_wall = (abs(rho(k)-R0)<tol || z(k)<tol || abs(z(k)-L0)<tol || \
	     abs(rho(k)-Rin)<tol);
  if !is_wall
    continue;			# for efficiency
  elseif !closed_tube && abs(rho(k)-R0)<tol && z(k)<=Dz_in
    ## Inlet
    er = x3d(k,[1 2]);
    er = er/l2(er);
    et = [-er(2) +er(1)];
    u = -u_rad_in * er + u_circunf_in * et;
    if inlets && compressible; fprintf(fid,"%d %d    %f\n",k,1,rho_in); endif
    fprintf(fid,"%d %d    %f\n",k,u_dof,u(1));
    fprintf(fid,"%d %d    %f\n",k,u_dof+1,u(2));
    fprintf(fid,"%d %d    %f\n",k,u_dof+2,0);
    n_in = n_in+1;
  elseif inlets && !closed_tube && abs(rho(k)-R0)<tol && z(k)>=L0-Dz_h
    fprintf(fid,"%d %d   %f\n",k,p_dof,p_h);
    fprintf(fid,"%d %d   %f\n",k,u_dof+2,0);
    n_h = n_h+1;
  elseif 0 && inlets && !closed_tube && z(k)<tol && rho(k)<=Rc
    fprintf(fid,"%d %d   %f\n",k,p_dof,p_c);
    n_c = n_c+1;
  elseif is_wall
    n_wall = n_wall+1;
    for l=1:3; fprintf(fid,"%d %d    %f\n",k,u_dof+l-1,0); endfor
  endif
endfor
fclose(fid);
printf(["Nodes at boundary %d, at inlet %d, at hot outlet %d," \
	" at cold outlet %d, at wall %d\n"],
       nn,n_in,n_h,n_c,n_wall);

## axisymmetric constraints
if axisymm
  fid = fopen("vtube.peri.tmp","w");
  c = cos(Dphi);
  s = sin(Dphi);
  for k=1:nn
    peri = k+nn;
    fprintf(fid,"%f  %d %d     %f  %d %d       %f  %d %d\n",
	    -1.,peri,u_dof  ,+c,k,u_dof,-s,k,u_dof+1);
    fprintf(fid,"%f  %d %d     %f  %d %d       %f  %d %d\n",
	    -1.,peri,u_dof+1,+s,k,u_dof,+c,k,u_dof+1);
    odof = [1 4 5];
    if !compressible, odof = [3 4]; endif
    for l=odof
      fprintf(fid,"%f  %d %d     %f  %d %d\n", \
	      -1.,peri,l,+1.,k,l);
    endfor
  endfor
  fclose(fid);
endif

Omega = u_circunf_in/R0;
nnod = rows(x3d);
if !compressible, p_in=0; endif
uini = ones(nnod,4)*diag([0,0,0,p_in]);
uini(:,1:2) = [-Omega*x3d(:,2),+Omega*x3d(:,1)];
if compressible
  uini = [rho_in*ones(nnod,1) uini];
endif
asave("vtube.ini.tmp",uini);
