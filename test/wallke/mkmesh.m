##__INSERT_LICENSE__
## $Id: mkmesh.m,v 1.7 2003/01/08 15:49:05 mstorti Exp $
startup
source("data.m.tmp");

hav=L/N;
w=zhomo([0 L 0 hav],N+1,2,[1 0 hratio 1 0 1]);
[xnod,icone]=pfcm2fem(w);
xnod = [xnod;
        xnod(1:2,:)];

icone=icone(:,[1 4 3 2]);
fid = fopen("wallke.peri.tmp","w");
for k=1:N+1
  for kd=1:5
    fprintf(fid,"%f %5d %5d %f %5d %5d\n",-1,2*k,kd,+1,2*k-1,kd);
  endfor
endfor
fclose(fid);

uini=[0. 0.666 0. 0.1 0.05];
uini=[uini(ones(2*(N+1),1),:);
      zeros(2,5)];

asave("wallke.nod.tmp",xnod);
asavecon("wallke.con.tmp",icone);
asave("wallke.ini.tmp",uini);

yp=0.05;
uini(:,2)=(yp+xnod(:,1))/(1+yp);;
asave("wallke.inic.tmp",uini);
