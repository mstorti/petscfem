## $Id: proc22d.m,v 1.5 2005/01/29 16:51:34 mstorti Exp $

source("data.m.tmp");

field = 1;

U = aload("gfabso2dn.some-rslt.tmp");
U(:,1)=[];
# gasdata.gamma = gamma;
# gasdata.pref = pref;
# gasdata.uref = uref;
# gasdata.rhoref = rhoref;
# gasdata.cref = cref;

scale = [rhoref cref cref pref];

if 1
  ## Add Mach
  c = sqrt(gamma*U(:,1)./U(:,4));
  Ma = U(:,2)./c;
  U = [U,Ma];
  scale = [scale,1];
endif

## Do not scale
scale = ones(size(scale));

ndof = columns(U);
Uref = mean(U);
Uref = zeros(1,ndof);
maxbound = max(U);
minbound = min(U);

m = columns(U);
Unorm = zeros(rows(U),m);
for k=1:m
  Unorm(:,k) = (U(:,k)-Uref(k))/scale(k);
endfor

rem(rows(U),Nx+1)==0 || error("not correct size");
nt = rows(U)/(Nx+1);
x = aload("gfabso2dn.nod.tmp");
x = x(1:Nx+1,longindx);

axis([min(x) max(x) min(min(Unorm)) max(max(Unorm))]);
m=3;

for k=1:m:nt
  title(sprintf("step %d",k));
  plot(x,Unorm((k-1)*(Nx+1)+(1:Nx+1),:))
  pause(0.1);
  ## pause;
endfor
