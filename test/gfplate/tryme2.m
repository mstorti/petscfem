## $Id: tryme2.m,v 1.1 2005/01/23 18:25:12 mstorti Exp $

nor = 2*rand(1,3)-1;

nor = nor/l2(nor);
[bid,indx]=min(nor);
z = zeros(1,3);
z(indx)=1;

t1 = pvec(nor,z);
t1 = t1/l2(t1);
t2 = pvec(nor,t1);

O = [nor;t1;t2];
