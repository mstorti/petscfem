## $Id: proc2.m,v 1.2 2005/01/21 03:14:30 mstorti Exp $

source("data.m.tmp");

U = aload("gfabso.some-rslt.tmp");

rem(rows(U),Nx+1)==0 || error("not correct size");

nt = rows(U)/(Nx+1);
u=reshape(U(:,3),Nx+1,nt);
nnod = rows(u);
u = u(1:nnod-2,:);

plot(u)

