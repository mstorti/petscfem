##__INSERT_LICENSE__
## $Id: proc.m,v 1.4 2003/01/08 15:49:04 mstorti Exp $
source("data.m.tmp");

U=aload("save.state.tmp");
xnod=aload("plate.nod.tmp");
icone=aload("plate.con.tmp");
icone=[icone(:,[1 2 3]);
	icone(:,[3 4 1])];
nnod = length(xnod);
y=xnod(1:Nx+1:nnod,2);
x=xnod(1:Nx+1,1);
u=reshape(U(:,1),Nx+1,Ny+1)';
v=reshape(U(:,2),Nx+1,Ny+1)';
p=reshape(U(:,3),Nx+1,Ny+1)';

if 0
  plot(y,u);
  pause
  plot(y,v);
  pause
  plot(y,p);
  pause
endif

W = curl(xnod,icone,U(:,1:2));
w = reshape(W,Nx+1,Ny+1)';
plot(y,w);
