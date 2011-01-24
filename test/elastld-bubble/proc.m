###key proc.m
### $Id: proc.m,v 1.2 2007/01/26 02:37:55 mstorti Exp $
source("data.m.tmp");
icone = aload("bubble.con.tmp");
x0 = aload("bubble.nod.tmp");
agraph = gplfem(icone);

kinc = 3;
k = 0;

xh = [];
while 1
  file = sprintf("./STEPS/bubble.state-%d.tmp",k);
  if !existfile(file); break; endif;
  printf("processing file %s\n",file);
  u = aload(file);
  x = x0+u;
  xh = [xh;-x(1,1)];
  gplot(agraph,x);
  axis([0,1.5,0,1.5]*1.0);
  pause(0.1);
  k += kinc;
endwhile
