## $Id: proc2.m,v 1.4 2005/01/21 17:10:47 mstorti Exp $

source("data.m.tmp");

field = 1;

U = aload("gfabso.some-rslt.tmp");
rem(rows(U),Nx)==0 || error("not correct size");
nt = rows(U)/Nx;
u=reshape(U(:,3),Nx,nt);
nnod = rows(u);
x = aload("gfabso.nod.tmp");
x = x([1:Nx-1,Nx+1],1);

gasdata.gamma = gamma;
Uri = primi2ri(u,gasdata);

axis([0 Lx min(min(u)) max(max(u))])
m=1;
for k=1:m:columns(u)
  title(sprintf("step %d",k));
  plot(x,u(:,k));
  pause(0.1);
endfor
