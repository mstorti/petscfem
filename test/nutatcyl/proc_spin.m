##__INSERT_LICENSE__
## $Id: proc_spin.m,v 1.3 2003/01/08 15:49:04 mstorti Exp $
f = aload("cylinder.all_force_spin.tmp");

nomega = floor(rows(f)/2);
f=f(1:2*nomega,:);
nuta = reshape(f(:,2),2,nomega)';
(all(nuta(:,1)==0) && all(nuta(:,2)==15)) || error("not correct nutation angles");

Omega = reshape(f(:,1),2,nomega)';
Omega = f(1+(0:nomega-1)'*2,1);
[Omega,indx] = sort(Omega);

## We multiply the  computed moments by the density of the fluid
## (silicone oil 950Kg/m3) and a conversin factor lbf-ft/Newton-m = 0.738
Mz_raw = reshape(f(:,8),2,nomega)'*950*0.738;
Mz_raw = Mz_raw(indx,:);
Mz = -(Mz_raw(:,2)-Mz_raw(:,1));

