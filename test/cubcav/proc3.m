## $Id: proc3.m,v 1.1 2003/11/28 03:13:13 mstorti Exp $
source("data.m.tmp");

u=aload("cubcav.state-plain.tmp");

x = aload("cubcav.nod.tmp");
x = x(1:N+1,1);

for k=1:3
  proc4 (u(:,k),N,x);
endfor
