###key proc.m
### $Id: proc.m,v 1.2 2007/01/26 02:37:55 mstorti Exp $
source("data.m.tmp");
icone = aload("gascont.con.tmp");
icone = icone(:,1:3);
agraph = gplfem(icone);

kinc = 4;
k = 1;

while 1
  file = sprintf("./STEPS/gascont-mmv.state-%d.tmp",k);
  if !existfile([file ".gz"]); break; endif;
  printf("processing file %s\n",file);
  xnod = aload(file);
  xnod = xnod(:,1:2);
  gplot(agraph,xnod);
  axis([0,1,0,1]);
  pause(0.1);
  k += kinc;
endwhile
