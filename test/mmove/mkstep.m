source("data.m.tmp");
w = zhomo([0 1 0 1],N+1,N+1);
[xnod,icone] = pfcm2fem(w);

## icone = icone(:,[1 4 3 2]);
icone=[icone(:,[1 3 2]);
       icone(:,[1 3 4])];

asave("step.nod.tmp",xnod);
asave("step.con.tmp",icone);

fid = fopen("step.fixa.tmp","w");
for k=1:N+1
  fprintf(fid,"%d %d %f\n",k,1,0.);
  fprintf(fid,"%d %d %f\n",k,2,0.);
  fprintf(fid,"%d %d %f\n",(N+1)*N+k,1,0.);
  fprintf(fid,"%d %d %f\n",(N+1)*N+k,2,0.);

  if k!=1 && k!=N+1
    fprintf(fid,"%d %d %f\n",(N+1)*(k-1)+1,1,0.);
    fprintf(fid,"%d %d %f\n",(N+1)*(k-1)+1,2,0.);
    fprintf(fid,"%d %d %f\n",(N+1)*k,1,0.);
    fprintf(fid,"%d %d %f\n",(N+1)*k,2,0.);
  endif

endfor
fclose(fid);
