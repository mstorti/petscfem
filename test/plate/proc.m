source("data.m.tmp");

U=aload("save.state.tmp");
u=reshape(U(:,1),Nx+1,Ny+1)';
v=reshape(U(:,2),Nx+1,Ny+1)';
p=reshape(U(:,3),Nx+1,Ny+1)';

