if 0
  u = aload("./comp_corner.dx-state.tmp");
  x = aload("./comp-corner-mesh-57k/comp_corner_Ma_5_low_Rey.nod.tmp");
endif

nnod = rows(x);
Ny = 406;
windx = (1:(Ny+1):nnod)';
windx(length(windx))=[];

Cp = (u(windx,4)-u(1,4))./(0.5*u(1,1)*u(1,2)^2);
xw = x(windx,1);
