## $Id: proc2.m,v 1.3 2005/01/21 15:25:25 mstorti Exp $

source("data.m.tmp");

if 1
  U = aload("gfabso.some-rslt.tmp");
  rem(rows(U),Nx+1)==0 || error("not correct size");
  nt = rows(U)/(Nx+1);
  u=reshape(U(:,3),Nx+1,nt);
  nnod = rows(u);
  u = u(1:nnod-2,:);
  x = aload("gfabso.nod.tmp");
  x = x(1:Nx-1,1);
endif

axis([0 Lx min(min(u)) max(max(u))])
m=1;
for k=1:m:columns(u)
  plot(x,u(:,k));
  pause(0.2);
endfor

