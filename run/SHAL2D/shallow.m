%% $Id: shallow.m,v 1.1.1.1 2000/12/28 12:54:43 mstorti Exp $
global gravity addvisc Lx Ly Nx Ny dx dy dHdx dHdy flag_1dx

if exist('restart') %&& restart
  answ=input('Continuar la corrida previa? (y/n) > ','s');
  if  ~strcmp(answ,'y') 
    return
  end
  nprev=size(res,1);
  res=[res;
       zeros(ntime,3)];
else

  gravity=1;
  Nx=80; Lx=8; 
  Ny=1; Ly=Ny*(Lx/Nx);
  ntime=100; Cou=0.7; 
  DH=.2; sigma=1; addvisc=.0;
  cfric=0.1;
  flag_1dx=1;
  NN=(Nx+1)*(Ny+1);

  % Posicion de la perturbacion
  Xc=Lx/2;

  % metodo de integracion temporal
  metodo='explicito';

  % fondo variable en el tiempo?
  fondovar=0;
  %omega=2*pi/(10*dt);
  omega=2;

  % salvar cada nsave, un maximo de nsavemax
  nsave=4; nsavemax=25;

  % salvar los valores en el vector IIsave
  %ky=1;
  %IIsave=(1:Nx+1)'+(ky-1)*(Nx+1);
  IIsave=(1:NN);

  % estados de un lado y del otro
  u0=[1 0.4 0];
  uL=u0; uL(1)=.3;
  
  cmax=sqrt(u0(2)^2+u0(3)^2)+sqrt(gravity*u0(1));

  x=(0:Nx)'/Nx*Lx;
  dx=(x(2)-x(1));
  y=(0:Ny)'/Ny*Ly;
  dy=(y(2)-y(1));

  H=zeros(Nx+1,Ny+1);;
  U0=[u0(1),u0(1)*u0(2),u0(1)*u0(3)];
  UL=[uL(1),uL(1)*uL(2),uL(1)*uL(3)];

  U = zeros(NN,3);
  for ky=1:Ny+1
    II=(1:Nx+1)'+(ky-1)*(Nx+1);
    r2=(x-Xc).^2+y(ky)^2;
    if flag_1dx, r2=(x-Xc).^2; endif % para hacerlo 1d
    H(:,ky)=H(:,ky)+DH*exp(-r2/(0.5)^2);
    U(II,:)=kron((U0+UL)/2,ones(size(x)))-kron((U0-UL)/2,tanh((x-Xc)/.01));
    % pone campana en la altura 
    %H=0*H;
    %U(II,1)=U(II,1)+DH*exp(-r2/(0.5)^2);
  end

  dHdx=[zeros(1,Ny+1);xcent(diff(H)/dx);zeros(1,Ny+1)];
  dHdx=dHdx(:);
  dHdy=[zeros(Nx+1,1),xcent(diff(H')/dx)',zeros(Nx+1,1)];
  dHdy=dHdy(:);

  if fondovar
    dHdx0=dHdx;
    dHdy0=dHdy;
  end

  dt=Cou*min([dx dy])/cmax
  R=zeros(size(U));

  % redondea a T=numero entero de pasos de tiempo
  if fondovar
    T=2*pi/omega;
    n=ceil(T/dt/nsave);
    dt=T/n/nsave;
    if T/dt/nsave>nsavemax
      warning('Atencion que el numero de pasos necesarios no entra en nsavemax!!')
    end
  end

  nprev=0;
  kt=0;
  res=[];
  hh=zeros(length(IIsave),nsavemax);
  uh=hh;
  vh=hh;

end

tic;
isave=0;
for it=1:ntime
  kt=kt+1;
  %fflush(stdout);
  if rem(kt-1,nsave)==0
    isave=isave+1;
    if isave>nsavemax
      isave=isave-nsavemax;
    end
    disp(sprintf('Save: kt %d, isave %d',kt, isave));
    hh(:,isave)=U(IIsave,1);
    uh(:,isave)=U(IIsave,2)./U(IIsave,1);
    vh(:,isave)=U(IIsave,3)./U(IIsave,1);
  end

  if fondovar
    dHdx=cos(omega*kt*dt)*dHdx0;
    dHdy=cos(omega*kt*dt)*dHdy0;
  end

  R=resshal(U);
  
  if strcmp('metodo','newton')
    Jaco=zeros(3*(Nx+1)*(Ny+1));
    jnod=0;  ldof=0;
    epsil=1e-6;
    for ky=1:Ny+1
      for kx=1:Nx+1
	jnod=jnod+1;
	for jdof=1:3
          ldof=ldof+1;
          UU=U;
          UU(jnod,jdof)=UU(jnod,jdof)+epsil;
          RR=resshal(UU);
          Jaco(:,ldof)=reshape((RR-R)'/epsil,3*NN,1);
	end
      end
    end
    dU=(eye(3*NN)+dt*Jaco)\(reshape(R',3*NN,1));
    UU=U-dt*reshape(dU,3,NN)';
  elseif strcmp(metodo,'explicito')
    UU=U-dt*R;
  else
    disp(sprintf('No existe metodo: %s',metodo));
  end

  resu=sqrt(sum((U-UU).^2));
  res(kt,:)=resu;
  disp(sprintf('Iter: %4d   Elaps: %5.1f  Res(h|u|v): %5.3e %5.3e %5.3e',kt,toc,resu));
  U=UU;
end
