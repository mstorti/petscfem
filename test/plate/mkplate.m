source("data.m.tmp");

w=zhomo([0 Ly 0 Lx],Ny+1,Nx+1,[1 0 yratio 1 0 1]);
[xnod,icone] = pfcm2fem(w);
xnod = xnod(:,[2 1]);

xnod(:,1) = xnod(:,1)-xini;
x = xnod(:,1);
y = xnod(:,2);

dy = 4*t/Lplate^2* (x>0. & x<=1.).*x.*(Lplate-x).*(1-y/Ly);
xnod = [x y+dy];

asave("plate.nod.tmp",xnod);
asave("plate.con.tmp",icone);

nnod = rows(xnod);
ini = [1 0 0];
ini = ini(ones(nnod,1),:);
asave("plate.ini.tmp",ini);

fid = fopen("plate.fixa.tmp","w");
for k=1:Nx+1
  ## On the plate
  fprintf(fid,"%d %d %f\n",k,2,0.);
  ## Plate extends only from 0 to 1
  xx = xnod(k,1);
  if xx>=0 && xx <=1
    fprintf(fid,"%d %d %f\n",k,1,0.);
  endif
endfor

V = aload("v.tmp");
for k=2:Nx+1
  ## On the external boundary
  fprintf(fid,"%d %d %f\n",Ny*(Nx+1)+k,2,V(k));
endfor

for k=2:Ny
  ## u=1, v=0 at the inlet
  fprintf(fid,"%d %d %f\n",(k-1)*(Nx+1)+1,1,1.);
  fprintf(fid,"%d %d %f\n",(k-1)*(Nx+1)+1,2,0.);
  ## v=0 at the outlet
  fprintf(fid,"%d %d %f\n",k*(Nx+1),2,0.);
endfor

fprintf(fid,"%d %d %f\n",Ny*(Nx+1)+1,1,1.);
fprintf(fid,"%d %d %f\n",Ny*(Nx+1)+1,2,V(1));
fprintf(fid,"%d %d %f\n",(Ny+1)*(Nx+1),2,0.);

for k=1:Ny+1
  ## p=0 at the outlet
  fprintf(fid,"%d %d %f\n",k*(Nx+1),3,0.);
endfor

fclose(fid);
