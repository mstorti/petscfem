## Number os elements in x,y
source("sine.data");

Nx=nx+1;
Ny=ny+1;

w=zhomo([0 Lx 0 Ly],nx+1,ny+1);
[xnod,icone]=pfcm2fem(w);
icone=icone(:,[1 4 3 2]);

yy=xnod(1:ny+1,2);
yyamp=sin(yy*pi/2);
fixa=[(1:ny+1)' ones(ny+1,1) yyamp];

fixa=[fixa;
      (Ny+1:Ny:Nx*Ny)' ones(nx,1) zeros(nx,1)];

fid=fopen("sine.fixa.tmp","w");
fprintf(fid,"%d %d %f\n",fixa');
fclose(fid);

asave("sine.nod.tmp",xnod);
asave("sine.con.tmp",icone);

## bcconv en todo el fondo y la tapa
bcconv = [nx*Ny+(1:Ny)';
          (nx*Ny:-Ny:Ny)'];

some = Nx*Ny;
asave2("sine.some.tmp",some);

nbc=rows(bcconv);
bcconv = [bcconv(1:nbc-1) bcconv(2:nbc)];
asave("sine.bcconv.tmp",bcconv);
