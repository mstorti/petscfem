##__INSERT_LICENSE__
## $Id: mkmesh.m,v 1.6 2003/01/08 15:49:04 mstorti Exp $
source("data.m.tmp");

hav=Ly/Ny;
w=zhomo([0 Lx 0 Ly],Nx+1,Ny+1,[1 0 1 1 hratio 1]);
[xnod,icone]=pfcm2fem(w);
icone=icone(:,[1 4 3 2]);
nnod2d=rows(xnod);
[icone,xnod]=extrude(xnod,icone,Nz,Lz/Nz);
nnod = rows(xnod);

fid = fopen("stabi.peri.tmp","w");
## Periodicity in the x direction
for l=1:Nz+1
  for k=1:Ny+1
    for kd=1:3
      fprintf(fid,"%f %d %d %f %d %d\n",-1,
              nnod2d*(l-1)+Nx*(Ny+1)+k,kd,+1,nnod2d*(l-1)+k,kd);
    endfor
  endfor
endfor

## Periodicity in the z direction
for l=1:Nx+1
  for k=1:Ny+1
    nodo=(l-1)*(Ny+1)+k;
    nodop=nodo+Nz*nnod2d;
    for kd=1:3
      fprintf(fid,"%f %d %d %f %d %d\n",-1,nodop,kd,+1,nodo,kd);
    endfor
  endfor
endfor
fclose(fid);
asave("stabi.nod.tmp",xnod);
asave("stabi.con.tmp",icone);

fid = fopen("stabi.fixaw.tmp","w");
for l=1:Nz
  for k=1:Nx
    nodo1=(l-1)*nnod2d+(k-1)*(Ny+1)+1;
    nodo2=nodo1+Ny;
    for kd=1:3
      fprintf(fid,"%d %d %f\n",nodo1,kd,0.);
      fprintf(fid,"%d %d %f\n",nodo2,kd,0.);
    endfor
  endfor
endfor
fclose(fid);

umax=1.;
cpert=0.5;
nmode=4;
urand=1;
uini = zeros(nnod,4);
uini(:,1) = umax*4*xnod(:,2).*(Ly-xnod(:,2))/Ly^2;
for ikx=1:nmode
  for iky=1:nmode
    for ikz=1:nmode
      coef = 3. / (ikx^2+iky^2+ikz^2);
      upertx = (2*rand-1)*cos(2*pi*ikx*xnod(:,1)/Lx) \
	  + (2*rand-1)*sin(2*pi*ikx*xnod(:,1)/Lx);
      uperty = (2*rand-1)*sin(iky*pi*xnod(:,2)/Ly);
      upertz = (2*rand-1)*cos(2*pi*ikz*xnod(:,3)/Lz) \
	  + (2*rand-1)*sin(2*pi*ikz*xnod(:,3)/Lz);
      uini(:,1) = uini(:,1) + cpert*coef*upertx.*uperty.*upertz;
    endfor
  endfor
endfor
clear upertx uperty upertz
asave("stabi.ini.tmp",uini);

some_nodes = [(1:Ny+1)';
	      (Ny/2+1:nnod2d:nnod)'];
asave("stabi.some.tmp",some_nodes);
