###key proc.m
### $Id: proc.m,v 1.1 2006/03/12 12:07:37 mstorti Exp $

source("data.m.tmp");
x0 = aload("elastld.nod.tmp");

kinc = 2;
kstart = 0;
kend = 0;
kframe = 0;
xscale = 1.0;
nnodmax = 1000;

nnod = rows(x0);
if nnod>nnodmax
  indx = randindx(nnodmax,nnod);
else
  indx = (1:nnod)';
endif

if kstart!=0
  printf("starting at k %d\n",kstart);
endif
k = kstart;
while 1
  if kend>0 && k>kend; break; endif
  file = sprintf("./STEPS/elastld.state-%d.tmp.gz",k);
  if !exist(file); break; endif
  printf("loading file %s\n",file);
  dx = aload(file);
  x = x0 + dx(:,1:3);
  plot(x(indx,1),x(indx,2),"o");
  axis([-0.6 0.6 -0.1 1.1]*Ly*1.);
  pause(0.4);
  k += kinc;
endwhile
