##__INSERT_LICENSE__
## $Id: smoke.m,v 1.3 2003/06/01 15:55:35 mstorti Exp $
source("data.m.tmp");

w=zhomo([0 Lx 0 Ly],N+1,M+1,[1 ratio 1 1 ratio 1]);
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

uini = 2*rand(N+1,M+1)-1;
uini = smsmooth(uini,n_smoth_steps,omega_smooth);
uini=vec(uini);
uimag = smsmooth(2*rand(N+1,M+1)-1,n_smoth_steps,omega_smooth);
uini = [uini vec(uimag)];
minu = min(min(uini));
maxu = max(max(uini));
uini=2*(uini-minu)/(maxu-minu)-1;
asave("smoke.ini.tmp",uini(:,1));
ue = pfnd2ele(xnod,icone,uini);

xele = pfnd2ele(xnod,icone,xnod(:,1));
yele = pfnd2ele(xnod,icone,xnod(:,2));
uele = [(yele-0.5).*xele.*(1-xele) -(xele-0.5).*yele.*(1-yele)];
fid = fopen("smoke.con.tmp","w");
for k=1:rows(icone)
  fprintf(fid,"%d %d %d    %f %f %f %f\n",icone(k,:),uele(k,:),ue(k,1),ue(k,2));
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
