## proc2

dir = "/u/jdelia/PETSC/CHANNEL3D";
source([dir "/data.m.tmp"]);
if 1,

  arch = [dir "/channel3d.nod.tmp"];
  x = load (arch) ;

  arch = [dir "/channel3d.some.tmp"];
  some = aload (arch);
  iy = some (1:Ny+1);
  iz = some;
  iz (1:Ny+1)=[];

  y = x (iy,2);
  z = x (iz,3);
  rem (Ny,2) == 0 || warning("Sould use Ny even!!");
  iz = [round(Ny/2)+1 Ny+1+(1:Nz)]';
  clear x

endif

nu = 0.000133;

if 0
nplane = (Nx+1)*(Nz+1);
nnod = nplane * (Ny+1);
uav= [];
uavh = zeros(Ny+1,1);
for k=0:30
  k
  u = aload([dir "/RUN/channel3d.state." int2str(k) ".tmp"]);
  for iy = 1:Ny+1;
    uavh(iy) = sum(u(iy+(0:(Ny+1):nnod-1),1))/nplane;
  endfor
  uav = [uav uavh];
endfor

endif
