source("sine.data");
u=aload("outvector0.sal");
nn = length(u);
N=(nx+1)*(ny+1);

nstep = round(nn/N);
if nstep*N != nn
  error("Not integer number of time steps");
endif
u=reshape(u,N,nstep);

Ix=(1:ny+1:N)';
xx=xnod(Ix,1);

Iy=(1:Ny+1)';
yy=xnod(Iy,2);
