##__INSERT_LICENSE__
## $Id: proc.m,v 1.4 2003/05/25 13:52:05 mstorti Exp $
source("data.m.tmp");
ghia;

load -force sqcav.ny.tmp
u=aload("sqcav.weak_form_0.tmp");
plot(u(ny,1),yh,u(ny,1),yh,'o',ug_1000,yg,'+')
