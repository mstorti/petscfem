## $Id: proc7.m,v 1.9 2005/03/02 12:08:24 mstorti Exp $

source("data.m.tmp");

if !exist("do_load") || do_load || !exist("U0");
  load -force ./gfshock3d.some-rslt.tmp
  U0 = gfshock3d;
  clear gfshock3d;
  ##  U0 = aload("gfshock3d.some-rslt.tmp");
  nod_some = unique(U0(:,1))';
  nsome = length(nod_some);
  all(U0(1:nsome,1)==nod_some) || error("bad some nodes");
  rem(rows(U0),nsome)==0 \
      || error("U0 not correct size.");
  nt = rows(U0)/nsome;
  U0(:,1)=[];
  if exist("discard") && discard>0
    printf("discarding %d steps\n",discard);
    U0 = U0((1:(nt-discard)*nsome),:);
    clear discard
  endif
endif
U=U0;
rem(rows(U),nsome)==0 || error("not correct size");
nt = rows(U)/nsome;

cin0 = sqrt(ga*pin0/rhoin0);
scale = [rhoin0 cin0 cin0 cin0 pin0];

ndof = columns(U);
Uref = zeros(1,ndof);
maxbound = max(U);
minbound = min(U);

if 0
  Uref = mean(U);
  scale=maxbound-minbound;
endif

m = columns(U);
Unorm = zeros(rows(U),m);
for k=1:m
  Unorm(:,k) = (U(:,k)-Uref(k))/scale(k);
endfor

x = aload("gfshock3d.nod.tmp");
x = x(nod_some,3);

Ue1 = U(nsome*(1:nt),:);
Ue0 = U(1+nsome*(0:nt-1),:);
rho = reshape(U(:,1),nsome,nt);

axis([min(x) max(x) min(min(Unorm)) max(max(Unorm))]);
m=1;

for k=1:m:nt
  title(sprintf("step %d",k));
  plot(x,Unorm((k-1)*(nsome)+(1:nsome),:))
  pause(0.1);
##  pause;
endfor
axis;
