##__INSERT_LICENSE__
## $Id: mkvtube.m,v 1.2 2003/01/17 19:03:00 mstorti Exp $

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

[x3d,ic3d] = extrude(xnod,icone,Nphi,1/Nphi);
z = x3d(:,1);
rho = x3d(:,2);
phi = 2*pi*x3d(:,3);

x3d = [rho.*cos(phi) rho.*sin(phi) z];

## Paste faces at phi=0, 2*pi
iren = (1:rows(x3d))';
iren(Nphi*nnod+(1:nnod)) = (1:nnod)';
for k=1:8
  ic3d(:,k) = iren(ic3d(:,k));
endfor
x3d(Nphi*nnod+(1:nnod),:) = [];

asave("vtube.nod.tmp",x3d);
asave("vtube.con.tmp",ic3d);

fid = fopen("vtube.fixa.tmp","w");
tol = 1e-7;
#wall = find((abs(rho-R0)<tol || z<tol || abs(z-L0)<tol)
done = 0;
for k=1:rows(x3d);
  if k/rows(x3d) > done+0.1;
    done = done+0.1;
    printf("%3d%% done\n",round(100*done));
  endif
  is_wall = (abs(rho(k)-R0)<tol || z(k)<tol || abs(z(k)-L0)<tol);
  if !is_wall
    continue;			# for efficiency
  elseif !closed_tube && abs(rho(k)-R0)<tol && z(k)<=Dz_in
    ## Inlet
    er = x3d(k,[1 2]);
    er = er/l2(er);
    et = [-er(2) +er(1)];
    u = -u_rad_in * er + u_circunf_in * et;
    fprintf(fid,"%d %d    %f\n",k,1,rho_in);
    fprintf(fid,"%d %d    %f\n",k,2,u(1));
    fprintf(fid,"%d %d    %f\n",k,3,u(2));
    fprintf(fid,"%d %d    %f\n",k,4,0);
  elseif !closed_tube && abs(rho(k)-R0)<tol && z(k)>=L0-Dz_h
    fprintf(fid,"%d %d   %f\n",k,5,p_h);
  elseif !closed_tube && z(k)<tol && rho(k)<=Rc
    fprintf(fid,"%d %d   %f\n",k,5,p_c);
  elseif is_wall
    fprintf(fid,"%d %d    %f\n",k,2,0);
    fprintf(fid,"%d %d    %f\n",k,3,0);
    fprintf(fid,"%d %d    %f\n",k,4,0);
  endif
endfor
fclose(fid);

Omega = u_circunf_in/R0;
nnod = rows(x3d);
uini = ones(nnod,5)*diag([rho_in,0,0,0,p_in]);
uini(:,2:3) = [-Omega*x3d(:,2),+Omega*x3d(:,1)];
asave("vtube.ini.tmp",uini);
