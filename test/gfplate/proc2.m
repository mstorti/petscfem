## $Id: proc2.m,v 1.7 2005/01/22 22:10:20 mstorti Exp $

source("data.m.tmp");

field = 1;

Uprimi = aload("gfabso.some-rslt.tmp");
Uprimi(:,1)=[];
gasdata.gamma = gamma;
gasdata.pref = pref;
gasdata.uref = uref;
gasdata.rhoref = rhoref;
gasdata.cref = cref;

primi = 0;
if primi
  U = Uprimi;
  field = 4;
  indx = [1 2 4];
  scale = [rhoref cref pref];
else
  Uri = primi2ri(Uprimi,gasdata);
  U = Uri;
  field = 1;
  indx = [1 2 3];
  scale = [cref cref 10.0];
endif

Uref = mean(U);
maxbound = max(U);
minbound = min(U);

Unorm = zeros(rows(U),3);
for k=1:3
  field = indx(k);
  Unorm(:,k) = (U(:,field)-Uref(field))/scale(k);
endfor

rem(rows(U),Nx+1)==0 || error("not correct size");
nt = rows(U)/(Nx+1);
x = aload("gfabso.nod.tmp");
x = x(1:Nx+1,1);

axis([0 Lx min(min(Unorm)) max(max(Unorm))])
m=5;

for k=1:m:nt
  title(sprintf("step %d",k));
  plot(x,Unorm((k-1)*(Nx+1)+(1:Nx+1),:))
  pause(0.1);
endfor
