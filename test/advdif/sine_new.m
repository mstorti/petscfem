## Number os elements in x,y
source("sine.data");
ndof=3;

Nx=nx+1;
Ny=ny+1;

w=zhomo([0 Lx 0 Ly],nx+1,ny+1);
[xnod,icone]=pfcm2fem(w);
icone=icone(:,[1 4 3 2]);

yy=xnod(1:ny+1,2);

fixa=[];
for j=1:3
  fixa=[fixa;
        (1:ny+1)' j*ones(ny+1,1) ones(ny+1,1)];
#    fixa=[fixa;
#          (1:ny+1)' j*ones(ny+1,2)];
endfor

fid=fopen("sine_new.fixa.tmp","w");
fprintf(fid,"%d %d %f\n",fixa');
fclose(fid);


nele=rows(icone);
if per_elem_prop
  xe=zeros(nele,1);
  for k=1:4
    xe=xe+xnod(icone(:,k),1);
  endfor
  xe=xe/4;
endif

asave("sine.nod.tmp",xnod);
#asave("sine.con.tmp",icone);
fid=fopen("sine.con.tmp","w");
for k=1:nele
  fprintf(fid,"%d %d %d %d",icone(k,1),icone(k,2),icone(k,3),icone(k,4));
  if per_elem_prop
#    fprintf(fid," %f %f ",xe(k),2*xe(k));
  endif
  fprintf(fid,"\n");
endfor
fclose(fid);

## bcconv en todo el fondo y la tapa
bcconv = [nx*Ny+(1:Ny)';
          (nx*Ny:-Ny:Ny)'];

some = Nx*Ny;
asave("sine.some.tmp",some);

nbc=rows(bcconv);
bcconv = [bcconv(1:nbc-1) bcconv(2:nbc)];
asave("sine.bcconv.tmp",bcconv);
