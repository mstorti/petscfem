source("data.m.tmp");

dteta = pi/100;
hav=L/N;
w=zhomo([L 0 0 hav],N+1,2,[1 0 hratio 1 0 1]);
[x2,ico2]=pfcm2fem(w);
x2 = [x2;
      x2(1:2,:)];

[icone,xnod]= extrude (x2,ico2,1,dteta);

Nlay = 2*N+4;
fid = fopen("wallke.peri.tmp","w");
fidf = fopen("wallke.fixa.tmp","w");
f="%f %5d %5d        %f %5d %5d\n";
for k=1:N+2
  for kd=1:6

    fprintf(fid,f,-1,       2*k,kd,+1,2*k-1,kd);
    fprintf(fid,f,-1,Nlay+2*k-1,kd,+1,2*k-1,kd);
    fprintf(fid,f,-1,  Nlay+2*k,kd,+1,2*k-1,kd);

    fprintf(fidf,"%5d %5d %f\n",       2*k,kd,0.);
    fprintf(fidf,"%5d %5d %f\n",Nlay+2*k-1,kd,0.);
    fprintf(fidf,"%5d %5d %f\n",  Nlay+2*k,kd,0.);

  endfor
endfor
fclose(fid);
fclose(fidf);

nnod=rows(xnod);
uini=[0. 0.666 0. 0. 0.1 0.05];
uini=[uini(ones(nnod,1),:)];

asave("wallke.nod.tmp",xnod);
asave("wallke.con.tmp",icone);
asave("wallke.ini.tmp",uini);
