source("data.m.tmp");
omega = 0.2;

## Initialization
if 0
  V = zeros(Nx+1,1);
  Vold = V;
  Vh = [];
endif

while 1
  V = (1-omega)*Vold+omega*V;
  asave("v.tmp",V);
  system("make run");

  Vold = V;
  Vh = [Vh V];
  plot(Vh);
  U=aload("save.state.tmp");
  u=reshape(U(:,1),Nx+1,Ny+1)';

  U=u(Ny+1,:)';
  U=[U;
     U(Nx:-1:1)];

  hx = Lx/Nx;
  VV = zeros(size(U));
  nnodx = rows(U);

  ## Loop over source panels
  for k=1:nnodx
    xrel = ((1:nnodx)'-k)*hx;
    ## Resolves singularity
    xrel(k) = 1;
    ixrel = 1./xrel;
    ixrel(k) = 0;
    VV = VV -(U(k)-1)*hx/pi*ixrel;
  endfor
  
  V=VV(1:Nx+1);
endwhile
