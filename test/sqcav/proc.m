##__INSERT_LICENSE__
## $Id: proc.m,v 1.5 2004/11/11 21:52:26 mstorti Exp $
source("data.m.tmp");
ghia;

load -force sqcav.ny.tmp
u=aload("sqcav.fractional_step_re1000.tmp");
if 1
  usome = aload("sqcav.some-state.tmp");
  nny = rows(ny);
  rem(rows(usome),nny)==0 || error("usome has incompatible size\n");
  ncol = rows(usome)/nny;
  usome = reshape(usome(:,2),nny,ncol);
  ncol = columns(usome);
  plot([u(ny,1),usome(:,1:5:ncol)],yh,u(ny,1),yh,'o',ug_1000,yg,'+')
else
  plot(u(ny,1),yh,u(ny,1),yh,'o',ug_1000,yg,'+')
endif
