function R=resshal(U)

  global Lx Ly Nx Ny dx dy gravity dHdx dHdy flag_1dx

  R=zeros(size(U));
  % Diferencia de flujos en la direccion X
  for ky=1:Ny+1
    II=(1:Nx+1)'+(ky-1)*(Nx+1);
    flux=upflux2(U([II(1);II;II(Nx+1)],:),[1 0]);
    R(II,:)=R(II,:)+(flux(2:Nx+2,:)-flux(1:Nx+1,:))/dx;
    R(II,2)=R(II,2)+gravity*U(II,1).*dHdx(II);
  end
  
  % Diferencia de flujos en la direccion Y
  if ~flag_1dx
    for kx=1:Nx+1
      II=kx+(Nx+1)*(0:Ny)';
      Uex0=U(II(2),:);
      Uex0(3)=-Uex0(3);
      UexL=U(II(Ny),:);
      UexL(3)=-UexL(3);
      flux=upflux2([Uex0;U(II,:);UexL],[0 1]);
      R(II,:)=R(II,:)+(flux(2:Ny+2,:)-flux(1:Ny+1,:))/dy;
      R(II,3)=R(II,3)+gravity*U(II,1).*dHdy(II);
    end
  end

  for ky=1:Ny+1
    II=1+(ky-1)*(Nx+1);
    R(II,:)=absorb2(R(II,:),U(II,:),[1 0]);
    II=Nx+1+(ky-1)*(Nx+1);
    R(II,:)=absorb2(R(II,:),U(II,:),[-1 0]);
  end
