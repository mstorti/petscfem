source("data.m.tmp");
ghia;

load -force sqcav.ny.tmp
u=aload("outvector0.out");
plot(u(ny,1),yh,u(ny,1),yh,'o',ug_1000,yg,'+')
