##__INSERT_LICENSE__
## $Id: smoke.m,v 1.2 2003/05/26 03:08:09 mstorti Exp $
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

uini = rand(N+1,M+1);
u2 = uini;
for j=1:n_smoth_steps
  uini(:,M+1) = uini(:,1);
  uini(N+1,:) = uini(1,:);

  for i=1:M+1;
    ip=modulo(i,N)+1;
    im=modulo(i-2,N)+1;
    u2(i,:) = omega_smooth*uini(i,:)+(1-omega_smooth)/2*(uini(ip,:)+uini(im,:));
  endfor
  uini = u2;

  for i=1:N+1;
    ip=modulo(i,M)+1;
    im=modulo(i-2,M)+1;
    u2(:,i) = omega_smooth*uini(:,i)+(1-omega_smooth)/2*(uini(:,ip)+uini(:,im));
  endfor
  uini = u2;

endfor
uini=vec(uini);
uini=(uini-min(uini))/(max(uini)-min(uini));
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
