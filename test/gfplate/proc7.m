## $Id: proc7.m,v 1.6 2005/02/26 15:59:56 mstorti Exp $

source("data.m.tmp");
nsome = Nx+1;

if !exist("do_load") || do_load || !exist("U0");
##  U0 = aload("gfabso2dn.some-rslt.tmp");
  U0 = aload("gfshock.some-rslt.tmp");
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
scale = [rhoin0 cin0 cin0 pin0];

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

x = aload("gfshock.nod.tmp");
x = x(1:nsome,1);

Ue1 = U(nsome*(1:nt),:);
Ue0 = U(1+nsome*(0:nt-1),:);
rho = reshape(U(:,1),nsome,nt);

axis([min(x) max(x) min(min(Unorm)) max(max(Unorm))]);
m=10;

for k=1:m:nt
  title(sprintf("step %d",k));
  plot(x,Unorm((k-1)*(nsome)+(1:nsome),:))
  pause(0.1);
##  pause;
endfor
axis;
