source("data.m.tmp");
w = zhomo([0 1 0 1],N+1,N+1,[4 0 1 1 0 1]);
[xnod,icone] = pfcm2fem(w);

  icone = icone(:,[1 4 3 2]);
if !strcmp(shape,"quad")
  icone=[icone(:,[2 3 1]);
	 icone(:,[4 1 3])];
endif

asave("step.nod.tmp",xnod);
asave("step.con.tmp",icone);

fid = fopen("step.fixa.tmp","w");
fidc = fopen("step.constrs.tmp","w");
for k=1:N+1
  y = xnod(k,2);
  if y<.3
    dd = 0.;
  elseif y>.7
    dd = disp;
  else
    dd = disp*(y-.3)/.4;
  endif
  dd = disp;
  fprintf(fid,"%d %d %f\n",k,1,dd);
  fprintf(fid,"%d %d %f\n",k,2,0.);
  fprintf(fid,"%d %d %f\n",(N+1)*N+k,1,0.);
  fprintf(fid,"%d %d %f\n",(N+1)*N+k,2,0.);

  if k!=1 && k!=N+1
    node1 = (N+1)*(k-1)+1;
    node2 = (N+1)*k;
    fprintf(fidc,"%f %d %d   %f %d %d\n",-1,node1,1,1,node2,1);
    fprintf(fidc,"%f %d %d   %f %d %d\n",-1,node1,2,1,node2,2);
  endif

endfor
fclose(fid);
fclose(fidc);

piecewtanh(slope,"piecewise.dat.tmp");
