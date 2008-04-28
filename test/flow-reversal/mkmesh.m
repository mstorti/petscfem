## __INSERT_LICENSE__
## $Id: mkstrip2d.m,v 1.1 2003/01/09 02:37:45 mstorti Exp $

source("data.m.tmp");

w=zhomo([0 Lx 0 Ly],Nx+1,Ny+1);

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

tmp = complement(right,[bot;top]);
pffixa("strip2d.wall.tmp",tmp,[1:2,4]);
pfperi("strip2d.peri.tmp",left,right,1:3);

icor = [right(1:Ny),right(2:Ny+1)];
icol = [(Ny+1:-1:2)',(Ny:-1:1)'];

asave("strip2d.flowrev-con.tmp",[icor;icol]);
