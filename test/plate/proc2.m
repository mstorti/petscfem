##__INSERT_LICENSE__
## $Id: proc2.m,v 1.3 2003/01/08 15:49:04 mstorti Exp $
source("data.m.tmp");
U=aload("save.state.tmp");
u=reshape(U(:,1),Nx+1,Ny+1)';

U=u(Ny+1,:)';
U=[U;
   U(Nx:-1:1)];

hx = Lx/Nx;
V = zeros(size(U));
nnodx = rows(U);

## Loop over source panels
for k=1:nnodx
  xrel = ((1:nnodx)'-k)*hx;
  ## Resolves singularity
  xrel(k) = 1;
  ixrel = 1./xrel;
  ixrel(k) = 0;
  V = V -(U(k)-1)*hx/pi*ixrel;
endfor

V=V(1:Nx+1);
asave("v.tmp",0.1*V);
