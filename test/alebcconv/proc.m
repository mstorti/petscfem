###key proc.m
### $Id: proc.m,v 1.12 2006/08/09 03:21:03 mstorti Exp $
source("data.m.tmp");

dire = "./STEPS";
x = aload("gfabso.nod.tmp");

x = x(1:Nx+1,1);

qindx = 1; ## density
kinc = 5;
if !exist("Uh")
  Uh = [];
  kstart = 1;
else 
  kstart = kinc*(columns(Uh)+1);
endif

k = kstart;
played = 0;
while 1
  file = sprintf("%s/gfabso.state-%d.tmp",dire,k);
  if !existfile(file); break; endif
  uu = aload(file);
  uu = uu(1:Nx+1,qindx);
  Uh = [Uh,uu];
  plot(x,uu);
  axis([0,Lx,0,1.5]);
  printf("step %d, time %f\n",k,k*Dt);
  pause(0.1);
  k = k + kinc;
  played = 1;
endwhile

if !played
  ncol = columns(Uh);
  for kk=1:ncol
    plot(x,Uh(:,kk));
    axis([0,Lx,0,1.5]);
    k = kk*kinc;
    printf("step %d, time %f\n",k,k*Dt);
    pause(0.1);
  endfor
endif

pause

nt = columns(Uh);
t = (0:nt-1)'*Dt;
plot(t,mean(Uh)');
