##__INSERT_LICENSE__
## $Id: proc4.m,v 1.2 2003/01/08 15:49:04 mstorti Exp $
x = aload("step3d.nod.tmp");
dx = aload("step3d.state.tmp");
x = x + dx;
## icone = aload("step3d.con-tet.tmp");
if 1
  checktetra(x,icone);
else
  asave("step3d.defo_nod.tmp",x);
endif
