###key proc.m
##__INSERT_LICENSE__
## $Id: proc2.m,v 1.2 2003/11/11 21:40:33 mstorti Exp $
source("data.m.tmp");

U=aload_any("pool.state_%d.tmp");

Nb = reshape(U(:,1),Nx+1,Ny+1);
CO = reshape(U(:,2),Nx+1,Ny+1);
CN = reshape(U(:,3),Nx+1,Ny+1);
COd = reshape(U(:,4),Nx+1,Ny+1);
CNd = reshape(U(:,5),Nx+1,Ny+1);
