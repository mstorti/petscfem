## $Id: proc2.m,v 1.14 2005/01/23 22:46:34 mstorti Exp $

source("data.m.tmp");

field = 1;

U = aload("gfabso.some-rslt.tmp");
U(:,1)=[];
gasdata.gamma = gamma;
gasdata.pref = pref;
gasdata.uref = uref;
gasdata.rhoref = rhoref;
gasdata.cref = cref;

U = U;
scale = [rhoref cref cref cref pref];

ndof = columns(U);
Uref = mean(U);

Unorm = zeros(size(U));
for k=1:ndof
  Unorm(:,k) = (U(:,k)-Uref(k))/scale(k);
endfor

rem(rows(U),Nx+1)==0 || error("not correct size");
nt = rows(U)/(Nx+1);
x = aload("gfabso.nod.tmp");
x = x(1:Nx+1,:)*nor';;

axis([min(x) max(x) \
      min(min(Unorm)) max(max(Unorm))]);
m = 5;

for k=1:m:nt
  k
  title(sprintf("step %d",k));
  plot(x,Unorm((k-1)*(Nx+1)+(1:Nx+1),:))
  pause(0.1);
#  pause;
endfor
