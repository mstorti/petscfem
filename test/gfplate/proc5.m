## $Id: proc5.m,v 1.10 2005/02/02 22:56:04 mstorti Exp $

source("data.m.tmp");

## line = "upstream";
## line = "downstream";
line = "axis";

U = aload("cylabso.some-rslt.tmp");
nline = Nphi+1;
nnod = nline*(Nr+1);
downstream = 1+nline*(0:Nr)';
upstream = nline*(1:Nr+1)';
xaxis = [nline*(Nr+1:-1:1),1+nline*(0:Nr)]';
skin = (1:Nphi+1);

U(:,1)=[];
gasdata.gamma = gamma;
gasdata.pref = pref;
gasdata.uref = uref;
gasdata.rhoref = rhoref;
gasdata.cref = cref;

ndof = columns(U);
U(:,ndof+1) = l2(U(:,2:ndof-1)) \
    ./sqrt(gamma*U(:,ndof)./U(:,1));

Uref = mean(U);

some = aload("cylabso.some-nodes.tmp");
nsome = length(some);

rem(rows(U),nsome)==0 || error("not correct size");
nt = rows(U)/nsome;
x = aload("cylabso.nod.tmp");

if strcmp(line,"upstream")
  plotindx = upstream;
  x = x(plotindx,1);
elseif strcmp(line,"downstream")
  plotindx = downstream;
  x = x(plotindx,1);
elseif strcmp(line,"axis")
  plotindx = xaxis;
  x = x(plotindx,1);
elseif strcmp(line,"skin")
  plotindx = skin;
  x = x(plotindx,1);
else
  error("bad case");
endif

if 0
  for k=1:columns(U)
    U(:,k) = U(:,k)-mean(U(:,k));
  endfor
endif

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
#  plot(x,U((k-1)*nsome+indx,:),x,U((k-1)*nsome+indx,:),'o')
  plot(x,U((k-1)*nsome+indx,:))
  pause(0.1)
#  pause;
endfor
