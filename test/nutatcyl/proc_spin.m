f = aload("cylinder.all_force_spin.tmp");

nomega = rows(f)/2;
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

