## $Id: proc2.m,v 1.1 2005/01/20 17:42:18 mstorti Exp $

source("data.m.tmp");

U = aload("gfabso.some-rslt.tmp");

rem(rows(U),Nx+1)==0 || error("not correct size");

nt = rows(U)/(Nx+1);
u=reshape(U(:,3),Nx+1,nt);

plot(u)

