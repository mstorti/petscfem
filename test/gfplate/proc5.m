## $Id: proc5.m,v 1.2 2005/01/24 18:36:10 mstorti Exp $

source("data.m.tmp");

field = 1;

U = aload("cylabso.some-rslt.tmp");
nline = Nphi+1;
nnod = nline*(Nr+1);
downstream = 1+nline*(0:Nr)';
upstream = nline*(1:Nr+1)';

plotindx = upstream;

U(:,1)=[];
gasdata.gamma = gamma;
gasdata.pref = pref;
gasdata.uref = uref;
gasdata.rhoref = rhoref;
gasdata.cref = cref;

ndof = columns(U);
Uref = mean(U);

some = aload("cylabso.some-nodes.tmp");
nsome = length(some);

rem(rows(U),nsome)==0 || error("not correct size");
nt = rows(U)/nsome;
x = aload("cylabso.nod.tmp");
x = x(plotindx,1);

for k=1:columns(U)
  U(:,k) = U(:,k)-mean(U(:,k));
endfor

axis([min(x) max(x) \
      min(min(U)) max(max(U))]);
m = 5;

indx = [];
for k=1:length(plotindx)
  j = find(some==plotindx(k));
  indx = [indx;j];
endfor

for k=1:m:nt
  title(sprintf("step %d",k));
  plot(x,U((k-1)*nsome+indx,:))
  pause(0.1);
#  pause;
endfor
