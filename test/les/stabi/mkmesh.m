##__INSERT_LICENSE__
## $Id: mkmesh.m,v 1.3 2003/01/08 15:49:04 mstorti Exp $
source("ini.data");

hav=Ly/N;
w=zhomo([0 3*hav 0 Ly],nlay+1,N+1,[1 0 1 1 hratio 1]);
[xnod,icone]=pfcm2fem(w);
icone=icone(:,[1 4 3 2]);
fid = fopen("stabi.peri.tmp","w");
for k=1:N+1
  for kd=1:3
    fprintf(fid,"%f %d %d %f %d %d\n",-1,nlay*(N+1)+k,kd,+1,k,kd);
  endfor
endfor
fclose(fid);

fid = fopen("stabi.fixaw.tmp","w");
for k=1:nlay
  fprintf(fid,"%d %d %f\n",(k-1)*(N+1)+1,1,0.);
  fprintf(fid,"%d %d %f\n",(k-1)*(N+1)+1,2,0.);
  fprintf(fid,"%d %d %f\n",k*(N+1),1,0.);
  fprintf(fid,"%d %d %f\n",k*(N+1),2,0.);
endfor
fclose(fid);

asave("stabi.nod.tmp",xnod);
asave("stabi.con.tmp",icone);
