###key proc.m
### $Id: proc.m,v 1.2 2007/01/26 02:37:55 mstorti Exp $
source("data.m.tmp");

kinc = 10;
k = 1;
x = (0:N)'/N;

icone = aload("./gascont.con.tmp");
agraph = gplfem(icone);

while 1
  file = sprintf("./STEPS/gascont.state-%d.tmp",k);
  if !existfile([file ".gz"]); break; endif;
  printf("processing file %s\n",file);
  u = aload(file);
  ## u = reshape(u,N+1,N+1)';

  file = sprintf("./STEPS/gascont-mmv.state-%d.tmp",k);
  if !existfile([file ".gz"]); break; endif;
  printf("processing file %s\n",file);
  xnod = aload(file);
  plotmesh3d([xnod(:,1:2),u],icone);
  pause;
  k += kinc;
endwhile
