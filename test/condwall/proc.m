## $Id: proc.m,v 1.1 2005/03/28 16:42:59 mstorti Exp $
source("data.m.tmp");

x = aload("condwall.nod.tmp");
u = aload("condwall.state.tmp");

w = (0:(Nx1+Nx2+1)')*(Ny+1)+1;

plot(x(w,1),u(w,3),x(w,1),u(w,3),'o');
