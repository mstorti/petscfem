## $Id: proc2.m,v 1.5 2005/01/21 18:23:08 mstorti Exp $

source("data.m.tmp");

field = 1;

Uprimi = aload("gfabso.some-rslt.tmp");
Uprimi(:,1)=[];
makegasdata.gamma = gamma;
Uri = primi2ri(Uprimi,gasdata);
U = Uri;

rem(rows(U),Nx)==0 || error("not correct size");
nt = rows(U)/Nx;
u=reshape(U(:,1),Nx,nt);
nnod = rows(u);
x = aload("gfabso.nod.tmp");
x = x([1:Nx-1,Nx+1],1);

axis([0 Lx min(min(u)) max(max(u))])
m=1;
for k=1:m:columns(u)
  title(sprintf("step %d",k));
  plot(x,u(:,k));
  pause(0.1);
endfor
