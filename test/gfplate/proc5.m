## $Id: proc5.m,v 1.1 2005/01/24 16:03:33 mstorti Exp $

source("data.m.tmp");

field = 1;

U = aload("cylabso.some-rslt.tmp");
nline = Nphi+1;
nnod = nline*(Nr+1);
upstream = 1+nline*(0:Nr)';
downstream = nline*(1:Nr)';

plotindx= downstream;

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
