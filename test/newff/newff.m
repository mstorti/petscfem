source("data.m");
ndof=3;

## Number of nodes/elements in x,y
Nx=nx+1;
Ny=ny+1;

w=zhomo([0 Lx 0 Ly],nx+1,ny+1);
[xnod,icone]=pfcm2fem(w);
icone=icone(:,[1 4 3 2]);

yy=xnod(1:ny+1,2);

## left to right peri
peri=[(ny+1)*nx+(1:ny+1)' (1:ny+1)'];
## upper to lower
nnod=(nx+1)*(ny+1);
peri=[peri;
      ny+(1:(ny+1):nnod)'  (1:(ny+1):nnod)'];
fid=fopen("newff.peri.tmp","w");
for j=1:ndof
  for k=1:rows(peri);
    fprintf(fid,"%d   %d   %d   %d   %d   %d\n", \
            -1,peri(k,1),j,1,peri(k,2),j);
  endfor
endfor
fclose(fid);

fixa=[];
for j=1:3
  fixa=[fixa;
        (1:ny+1)' j*ones(ny+1,1) ones(ny+1,1)];
endfor

fid=fopen("newff.fixa.tmp","w");
fprintf(fid,"%d %d %f\n",fixa');
fclose(fid);

nele=rows(icone);
if per_elem_prop
  xe=zeros(nele,2);
  for k=1:4
    xe=xe+xnod(icone(:,k),:);
  endfor
  xe=xe/4;
endif

##sour=cos(kx*xe(:,1)*2*pi/Lx).*cos(ky*xe(:,2)*2*pi/Ly);
kwave=[kx*2*pi/Lx ky*2*pi/Ly];
phase=kwave(1)*xe(:,1)+kwave(2)*xe(:,2);
phase=exp(i*phase);
sour=real(phase);

asave("newff.nod.tmp",xnod);
#asave("newff.con.tmp",icone);
fid=fopen("newff.con.tmp","w");
for k=1:nele
  fprintf(fid,"%d %d %d %d",icone(k,1),icone(k,2),icone(k,3),icone(k,4));
  if per_elem_prop
    fprintf(fid," %f %f ",sour(k));
  endif
  fprintf(fid,"\n");
endfor
fclose(fid);

## bcconv en todo el fondo y la tapa
bcconv = [nx*Ny+(1:Ny)';
          (nx*Ny:-Ny:Ny)'];

some = Nx*Ny;
asave("newff.some.tmp",some);

nbc=rows(bcconv);
bcconv = [bcconv(1:nbc-1) bcconv(2:nbc)];
asave("newff.bcconv.tmp",bcconv);

##--<*>---//---<*>---//---<*>---//---<*>---//---<*>---//
##--<*>---//--- ANALYTICAL SOLUTION -<*>---//---<*>---//
##--<*>---//---<*>---//---<*>---//---<*>---//---<*>---// 
uu=[u 0];
DD=[D 0; 0 D];

beta = R+i*kwave*uu'+kwave*DD*kwave';
phase=kwave(1)*xnod(:,1)+kwave(2)*xnod(:,2);
phase=exp(i*phase);
uana = s/beta*(1-exp(-beta*nstep*Dt));
Uana = real(uana*phase);
#Uana=reshape(Uana,17,17)';
asave("uanalyt.tmp",Uana(:,[1 1 1]));
