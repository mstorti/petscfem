## $Id: interp4.m,v 1.1 2005/02/25 20:04:33 mstorti Exp $

if 1
  data = "/u/mstorti/PETSC/tatuus/data/";
  xnod1 = aload([data "static_p_blade.nod"]);
  p1 = aload([data "static_p_blade.p"]);
  ico1 = aload([data "blade.con"]);
  xnod2 = aload([data "patran.nod"]);
endif

n2 = 70;				# node to check
dx = xnod1;
for k=1:3
  dx(:,k) = dx(:,k) - xnod2(n2,k);
endfor
dx = l2(dx);
[dmin,indx] = min(dx);
printf("node2 %d, dmin %f, u %f\n",n2-1,dmin,p1(indx));
