## $Id: proc3.m,v 1.2 2003/12/08 13:07:08 mstorti Exp $
source("data.m.tmp");

## u=aload("cubcav.state-plain.tmp");
u=aload("cubcav.state.plain_bupl0.tmp");

x = aload("cubcav.nod.tmp");
x = x(1:N+1,1);

for k=1:3
  proc4 (u(:,k),N,x);
endfor
