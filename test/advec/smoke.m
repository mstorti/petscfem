##__INSERT_LICENSE__
## $Id: smoke.m,v 1.1 2003/01/15 23:50:33 mstorti Exp $
source("data.m.tmp");

w=zhomo([0 Lx 0 Ly],N+1,M+1);
[xnod,icone]=pfcm2fem(w);
if !exist("noise")
  noise = 0.;
endif
hx = Lx/N;
hy = Ly/M;
nnod = rows(xnod);

icone=[icone(:,[1 4 3]);
       icone(:,[3 2 1])];
	
asave("smoke.nod.tmp",xnod);
asave("smoke.con.tmp",icone);

uini = 2*rand(nnod,1)-1;
asave("smoke.ini.tmp",uini);

xele = pfnd2ele(xnod,icone,xnod(:,1));
yele = pfnd2ele(xnod,icone,xnod(:,2));
uele = [(yele-0.5).*xele.*(1-xele) -(xele-0.5).*yele.*(1-yele)];
fid = fopen("smoke.con.tmp","w");
for k=1:rows(icone)
  fprintf(fid,"%d %d %d    %f %f\n",icone(k,:),uele(k,:));
endfor
fclose(fid);

fid = fopen("smoke.fixa.tmp","w");
for k=1:M+1
  node = k;
  fprintf(fid,"%d %d %f\n",node,1,0);
  fprintf(fid,"%d %d %f\n",node+N*(M+1),1,0);
endfor
for k=2:N
  node = (k-1)*(M+1)+1;
  fprintf(fid,"%d %d %f\n",node,1,0);
  node = (k-1)*(M+1)+M+1;
  fprintf(fid,"%d %d %f\n",node,1,0);
endfor
fclose(fid);
