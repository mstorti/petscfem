##__INSERT_LICENSE__
## $Id: proc.m,v 1.3 2003/01/08 15:49:05 mstorti Exp $
source("data.m.tmp");
ghia;

load -force sqcav.ny.tmp
u=aload("outvector0.out");
plot(u(ny,1),yh,u(ny,1),yh,'o',ug_1000,yg,'+')
