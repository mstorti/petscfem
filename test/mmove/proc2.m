##__INSERT_LICENSE__
## $Id: proc2.m,v 1.4 2003/01/08 15:49:04 mstorti Exp $
x = aload("step3d.nod.tmp");
dx = aload("step3d.state.tmp");
x = x + dx;
if 0
  asave("step3d.defo_nod.tmp",x);
else
  icone = aload("step3d.con-tet.tmp");
  [minv,maxv] = checktetra(x,icone);
  printf("All tetra volumes > 0 OK ? %d\n",minv>0);
  printf("Min/max vol %f/%f\n",minv,maxv);
  ratio = maxv/minv;
  max_ratio = 25;
  printf("Max/min ratio < max_ratio OK ? %d, ratio %f, max_ratio %f\n", \
	 ratio<max_ratio,ratio,max_ratio);
endif
