###key proc.m
##__INSERT_LICENSE__
## $Id: proc2.m,v 1.1 2003/11/11 15:41:13 mstorti Exp $
source("data.m.tmp");

U=aload("pool.state.tmp");

Nb = reshape(U(:,1),Nx+1,Ny+1);
CO = reshape(U(:,2),Nx+1,Ny+1);
CN = reshape(U(:,3),Nx+1,Ny+1);
COd = reshape(U(:,4),Nx+1,Ny+1);
CNd = reshape(U(:,5),Nx+1,Ny+1);
