source("ini.data");

hav=R/N;
w=zhomo([0 3*hav 0 R],nlay+1,N+1,[1 0 1 1 hratio 1]);
[xnod,icone]=pfcm2fem(w);
icone=icone(:,[1 4 3 2]);
nnod2d=rows(xnod);
[icone,xnod]=extrude(xnod,icone,3,hav);

fid = fopen("stabi.peri.tmp","w");
for l=1:nlay+1
  for k=1:N+1
    for kd=1:3
      fprintf(fid,"%f %d %d %f %d %d\n",-1,
              nnod2d*(l-1)+nlay*(N+1)+k,kd,+1,nnod2d*(l-1)+k,kd);
    endfor
  endfor
endfor

for l=1:nlay+1
  for k=1:N+1
    nodo=(l-1)*(N+1)+k;
    nodop=nodo+nlay*nnod2d;
    for kd=1:3
      fprintf(fid,"%f %d %d %f %d %d\n",-1,nodop,kd,+1,nodo,kd);
    endfor
  endfor
endfor
fclose(fid);
asave("stabi.nod.tmp",xnod);
asave("stabi.con.tmp",icone);

fid = fopen("stabi.fixaw.tmp","w");
for l=1:nlay
  for k=1:nlay
    nodo1=(l-1)*nnod2d+(k-1)*(N+1)+1;
    nodo2=nodo1+N;
    fprintf(fid,"%d %d %f\n",nodo1,1,0.);
    fprintf(fid,"%d %d %f\n",nodo1,2,0.);
    fprintf(fid,"%d %d %f\n",nodo2,1,0.);
    fprintf(fid,"%d %d %f\n",nodo2,2,0.);
  endfor
endfor
fclose(fid);
