##__INSERT_LICENSE__
## $Id: proc3.m,v 1.2 2003/01/08 15:49:04 mstorti Exp $
x0=aload("step3d.nod.tmp");
dx=aload("step3d.state.tmp");
x=x0+dx;

n=round(rows(x)^(1/3));

plot(X,Y)
