source("data.m.tmp");

u=aload("aquist.rot0.tmp");
u=reshape(u,(Nx+2)*(Ny+1),nstep);
uu=u((Nx+1)*(Ny+1)+(1:Ny+1)',50);

load -force u_stream.ref

erro = merr(uu-uref);
tol = 1e-3;
printf("Test OK? > %d, [error=%g, tol=%g]\n",erro<tol,erro,tol);
