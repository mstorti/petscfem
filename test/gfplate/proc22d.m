## $Id: proc22d.m,v 1.1 2005/01/24 22:01:13 mstorti Exp $

source("data.m.tmp");

field = 1;

Uprimi = aload("gfabso2d.some-rslt.tmp");
Uprimi(:,1)=[];
gasdata.gamma = gamma;
gasdata.pref = pref;
gasdata.uref = uref;
gasdata.rhoref = rhoref;
gasdata.cref = cref;

primi = 1;
if primi
  U = Uprimi;
  field = 4;
  indx = [1 2 3 4];
  scale = [rhoref cref cref pref];
else
  error("not more supported!!");
  Uri = primi2ri(Uprimi,gasdata);
  U = Uri;
  field = 1;
  indx = [1 2 3];
  scale = [cref cref 10.0];
endif

ndof = columns(U);
Uref = mean(U);
## Uref = zeros(1,ndof);
maxbound = max(U);
minbound = min(U);

m = length(indx);
Unorm = zeros(rows(U),m);
for k=1:m
  field = indx(k);
  Unorm(:,k) = (U(:,field)-Uref(field))/scale(k);
endfor

rem(rows(U),Nx+1)==0 || error("not correct size");
nt = rows(U)/(Nx+1);
x = aload("gfabso2d.nod.tmp");
x = x(1:Nx+1,longindx);

axis([min(x) max(x) min(min(Unorm)) max(max(Unorm))]);
m=5;

for k=1:m:nt
  title(sprintf("step %d",k));
  plot(x,Unorm((k-1)*(Nx+1)+(1:Nx+1),:))
  pause(0.1);
#  pause;
endfor
