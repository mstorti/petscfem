tcini

hav=Ly/N;
w=zhomo([0 Ly 0 hav],2,N+1,[1 0 10 1 0 1]);
[xnod,icone]=pfcm2fem(w);
#icone=icone(:,[1 4 3 2]);
xnod=xnod(:,[2 1]);
xnod=[xnod zeros(rows(xnod),1)];
fid = fopen("turbchanw.peri.tmp","w");
for k=1:N+1
  for kd=1:5
    fprintf(fid,"%f %d %d %f %d %d\n",-1,2*k,kd,+1,2*k-1,kd);
  endfor
endfor
fclose(fid);

asave("turbchanw.nod.tmp",xnod);
asave("turbchanw.con.tmp",icone);
