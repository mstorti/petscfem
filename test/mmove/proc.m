##__INSERT_LICENSE__
## $Id: proc.m,v 1.8 2003/12/08 01:02:29 mstorti Exp $
x0 = aload("step.nod.tmp");
dx = aload("step.state.tmp");
icone = aload("step.con.tmp");
x = x0+dx;
if 0
  gplfem(x,icone,"malla.gpl");
else
  area = l2(checktri(x,icone));
  mina = min(area);
  maxa = max(area);
endif

max_ratio = 500;
printf("Mesh OK (all areas >0) ? %d\nMin area %f, max area %f\n",mina>0,mina,maxa);
printf("Vol ratio < max allowed OK ? %d\nmax/min %f, max allowed\n",maxa/mina<max_ratio,max_ratio);
