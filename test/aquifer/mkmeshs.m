source("data.m.tmp");

x=[(0:Nx)'/Nx*Lx zeros(Nx+1,1)];

if 0 # || circular_path
  phi = x(:,1)/Lx*pi;
  x = Lx/pi*[cos(phi) sin(phi) 0*phi];
endif

ico = [(1:Nx)' (2:Nx+1)'];

ini=ones(rows(x),1);
asave("stream.nod.tmp",x);
asave("stream.con.tmp",ico);
asave("stream.ini.tmp",ini);
