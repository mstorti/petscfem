x0 = aload("step.nod.tmp");
dx = aload("step.state.tmp");
icone = aload("step.con.tmp");
x = x0+dx;
if 0
  gplfem(x,icone,"malla.gpl");
else
  [mina,maxa] = checktri(x,icone);
endif

max_ratio = 500;
printf("Mesh OK (all areas >0) ? %d\nMin area %f, max area %f\n",mina>0,mina,maxa);
printf("Vol ratio < max allowed OK ? %d\nmax/min %f, max allowed\n",maxa/mina<max_ratio,max_ratio);
