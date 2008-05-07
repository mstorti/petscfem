## __INSERT_LICENSE__
## $Id: mkstrip2d.m,v 1.1 2003/01/09 02:37:45 mstorti Exp $

source("data.m.tmp");

w=zhomo([0 h*Nx 0 L],Nx+1,Ny+1,[1 0 1 1 xratio 1]);

[x2,i2]=pfcm2fem(w);
i2 = i2(:,[1 4 3 2]);
asave("strip2d.nod.tmp",x2);
asave("strip2d.con.tmp",i2);

nnod = rows(x2);		# number of real nodes
nelem = rows(i2);
NN = Ny+1;			# number of nodes in one dimension

bot = (0:Nx)'*(Ny+1)+1;
top = bot+Ny;

left = (1:Ny+1)';
right = Nx*(Ny+1)+(1:Ny+1)';

pffixa("strip2d.wall.tmp",[bot;top],1:2);
pfperi("strip2d.peri.tmp",left,right,1:3);

icobot = [bot(1:Nx),bot(2:Nx+1),ones(Nx,4),zeros(Nx,3)];
icotop = [top(2:Nx+1),top(1:Nx),ones(Nx,4),zeros(Nx,3)];
asave("strip2d.visc-force-con.tmp",[icobot;icotop]);
