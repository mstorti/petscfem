%% $Id: shallowo.m,v 1.1 2000/12/28 12:54:43 mstorti Exp $
global gravity addvisc

if exist('restart') && restart
  answ=input('Continuar la corrida previa? (y/n) > ','s');
  if  ~strcmp(answ,'y') 
    return
  end
  nprev=size(res,1);
  res=[res;
       zeros(ntime,3)];
else

  gravity=1;
  Nx=20; Lx=2; 
  Ny=1; Ly=Ny*(Lx/Nx);
  ntime=20; Cou=0.7; 
  DH=.2; sigma=1; addvisc=.2;
  cfric=0.1;
  flag_1dx=1;

  % salvar cada nsave, un maximo de nsavemax
  nsave=1; nsavemax=50;

  % salvar los valores en el vector IIsave
  ky=round((Ny+1)/2);
  IIsave=(1:Nx+1)'+(ky-1)*(Nx+1);

  % estados de un lado y del otro
  u0=[1 0.4 0];
  ##  uL=supercr(u0(1:2));
  uL=[0.5 0.4 0];
  
  NN=(Nx+1)*(Ny+1);
  cmax=sqrt(u0(2)^2+u0(3)^2)+sqrt(gravity*u0(1));

  x=(0:Nx)'/Nx*Lx;
  dx=(x(2)-x(1));
  y=(0:Ny)'/Ny*Ly;
  dy=(y(2)-y(1));

  H=zeros(Nx+1,Ny+1);;
  U0=[u0(1),u0(1)*u0(2),u0(1)*u0(3)];
  UL=[uL(1),uL(1)*uL(2),uL(1)*uL(3)];
  for ky=1:Ny+1
    II=(1:Nx+1)'+(ky-1)*(Nx+1);
    r2=(x-Lx/2).^2+y(ky)^2;
    r2=(x-Lx/2).^2; % para hacerlo 1d
    H(:,ky)=H(:,ky)+DH*exp(-r2/(0.5)^2);
    U(II,:)=kron((U0+UL)/2,ones(size(x)))-kron((U0-UL)/2,tanh((x-Lx/2)/.6));
  end

  dHdx=[zeros(1,Ny+1);xcent(diff(H)/dx);zeros(1,Ny+1)];
  dHdx=dHdx(:);
  dHdy=[zeros(Nx+1,1),xcent(diff(H')/dx)',zeros(Nx+1,1)];
  dHdy=dHdy(:);

  dt=Cou*min([dx dy])/cmax;
  R=zeros(size(U));

  nprev=0;
  kt=0;
  res=[];
  hh=zeros(length(II),nsavemax);
  uh=hh;
  vh=hh;
end

tic;
isave=0;
for it=1:ntime
  kt=kt+1;
  fflush(stdout);
  if rem(kt-1,nsave)==0
    isave=isave+1;
    if isave>nsavemax
      isave=isave-nsavemax;
    end
    printf("Save: kt %d, isave %d\n",kt, isave);
    hh(:,isave)=U(IIsave,1);
    uh(:,isave)=U(IIsave,2)./U(IIsave,1);
    vh(:,isave)=U(IIsave,3)./U(IIsave,1);
  end

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

  UU=U-dt*R;
  return
  for ky=1:Ny+1
    II=1+(ky-1)*(Nx+1);
    UU(II,:)=absorb(UU(II,:),U(II,:),[1 0]);
    II=Nx+1+(ky-1)*(Nx+1);
    UU(II,:)=absorb(UU(II,:),U(II,:),[-1 0]);
  end
  
  resu=sqrt(sum((U-UU).^2));
  res(nprev+kt,:)=resu;
  printf("Iter: %4d   Elaps: %5.1f  Res(h|u|v): %5.3e %5.3e %5.3e\n",kt,toc,resu);
  U=UU;
end
