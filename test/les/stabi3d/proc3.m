## proc2
nu = 0.000133;

start = 200;
last = columns(u);
nsteps = (last-start+1);

nplane = (Nx+1)*(Nz+1);
nnod = nplane * (Ny+1);
uav= [];
uavh = zeros(Ny+1,1);
for k=0:17
  k
  u = aload([dir "/RUN/channel3d.state." int2str(k) ".tmp"]);
  for iy = 1:Ny+1;
    uavh(iy) = sum(u(iy+(0:(Ny+1):nnod-1),1))/nplane;
  endfor
  uav = [uav uavh];
endfor

