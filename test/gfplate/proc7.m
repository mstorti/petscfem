## $Id: proc7.m,v 1.2 2005/02/25 21:56:39 mstorti Exp $

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
scale = [1 1 1 1];

if 0
  ## Add Mach
  c = sqrt(ga*U(:,1)./U(:,4));
  Ma = U(:,2)./c;
  U = [U,Ma];
  scale = [scale,1];
  U = U(:,[1,2,4,5]);
  scale = scale([1,2,4,5]);
  U = [U,U(:,3)./U(:,1)];
  scale = [scale,1];
endif

## Do not scale
scale = ones(size(scale));

ndof = columns(U);
Uref = zeros(1,ndof);
maxbound = max(U);
minbound = min(U);

if 1
  Uref = mean(U);
  scale=maxbound-minbound;
endif
scale(3)=1;

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
m=1;

for k=1:m:nt
  title(sprintf("step %d",k));
  plot(x,Unorm((k-1)*(nsome)+(1:nsome),:))
  pause(0.1);
##  pause;
endfor
axis;
