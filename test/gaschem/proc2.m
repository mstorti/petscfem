###key proc.m
##__INSERT_LICENSE__
## $Id: proc2.m,v 1.8 2003/11/13 13:46:30 mstorti Exp $
source("data.m.tmp");

U=aload_any("pool.state_%d.tmp",0);
## U = aload("pool.state.tmp");

Nb = reshape(U(:,1),Nx+1,Ny+1)*Nb_scale;
CO = reshape(U(:,2),Nx+1,Ny+1);
CN = reshape(U(:,3),Nx+1,Ny+1);
COd = reshape(U(:,4),Nx+1,Ny+1);
CNd = reshape(U(:,5),Nx+1,Ny+1);

if !exist("xO_in")
  disp("doesn't exist xO_in, set to 0.2\n");
  xO_in = 0.2;
endif

## Compute saturation
xnod = aload("pool.nod.tmp");
y = xnod(:,2);
pg = patm+rho_liq*gravity*(max(y)-y);
pg = reshape(pg,Nx+1,Ny+1);

satO = COd./(KO*patm*xO_in);
