###key proc.m
##__INSERT_LICENSE__
## $Id: proc.m,v 1.2 2003/11/11 02:13:00 mstorti Exp $
source("data.m.tmp");

u=aload("outvector.out");

U = u(2*(1:Ny+1)-1,:);
x = aload("pool.nod.tmp");
y = x(2*(1:Ny+1)-1,2);
