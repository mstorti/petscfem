###key proc.m
##__INSERT_LICENSE__
## $Id: proc.m,v 1.3 2003/11/11 12:47:00 mstorti Exp $
source("data.m.tmp");

u=aload("outvector.out");

U = u(2*(1:Ny+1)-1,:);
x = aload("pool.nod.tmp");
y = x(2*(1:Ny+1)-1,2);

Uscale = sqrt(max(abs(U)).*min(abs(U)))';
Uscale = 10.^(round(log10(Uscale)));

U = rightscal(U,1./Uscale);
