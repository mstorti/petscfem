###key proc.m
##__INSERT_LICENSE__
## $Id: proc2.m,v 1.3 2003/11/12 12:47:02 mstorti Exp $
source("data.m.tmp");



U=aload_any("pool.state_%d.tmp",0);

Nb = reshape(U(:,1),Nx+1,Ny+1);
CO = reshape(U(:,2),Nx+1,Ny+1);
CN = reshape(U(:,3),Nx+1,Ny+1);
COd = reshape(U(:,4),Nx+1,Ny+1);
CNd = reshape(U(:,5),Nx+1,Ny+1);
