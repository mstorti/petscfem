### key proc4.m
### $Id: proc6.m,v 1.1 2005/02/03 21:16:34 mstorti Exp $

1;

function T = myresha(file,Rgas,Nphi,Nr)
  U=aload(file);
  T = U(:,4)./U(:,1)/Rgas;
  nnod = (Nphi+1)*(Nr+1);
  T = reshape(T(1:nnod),Nphi+1,Nr+1);
endfunction 

source("data.m.tmp");
xnod=aload("cylabso.nod.tmp");

nnod = (Nphi+1)*(Nr+1);
nfica = 2*(Nphi+1);

ficw = nnod+nfica+(1:Nphi)';
nficw = length(ficw);

ds = R*2*pi/Nphi;
q = U(ficw,1)/ds;
phi = (0:Nphi-1)'/Nphi*2*pi;

## By difference quotients
TT2 = myresha("STEPS/cylabso.state_18.tmp",Rgas,Nphi,Nr);
TT1 = myresha("STEPS/cylabso.state_17.tmp",Rgas,Nphi,Nr);
T = TT1;

T1 = T(1:Nphi,1);
T2 = T(1:Nphi,2);
T3 = T(1:Nphi,3);
dr = xnod(Nphi+2,1)-xnod(1,1);

rv = xnod(1+(0:Nr)*(Nphi+1),1);

kond = 1e-3;
q1 = -kond*(T2-T1)/dr;
r = xnod(1+(Nphi+1)*(0:2),1);
coef = cloud(r,1,2);
dTdr2 = (coef*[T1,T2,T3]')';
q2 = -kond*dTdr2;

if 0
  gamma = 1.4;
  Cp = Rgas/(gamma-1);
  temp_term = (TT2(1:Nphi,1)-TT1(1:Nphi,1))/(rho*Cp);
endif

plot(phi,[q,q1,q2]);
