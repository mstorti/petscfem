##__INSERT_LICENSE__
## $Id: mkmeshs.m,v 1.5 2003/01/08 15:49:03 mstorti Exp $
source("data.m.tmp");

x=[(0:Nx)'/Nx*Lx zeros(Nx+1,1)];
s = x(:,1);                     # arc length

if circular_path
  phi = x(:,1)/Lx*pi;
  x = Lx/pi*[cos(phi) sin(phi)];
endif

x=[x -slope*s];
ico = [(1:Nx)' (2:Nx+1)'];

ini=ones(rows(x),1);
asave("stream.nod.tmp",x);
asave("stream.con.tmp",ico);
asave("stream.ini.tmp",ini);
