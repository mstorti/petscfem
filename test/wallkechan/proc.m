##__INSERT_LICENSE__
## $Id: proc.m,v 1.11 2003/01/08 15:49:05 mstorti Exp $
source("data.m.tmp");
NN=(N+1);
NNy = Ny+1;
nnod = NN*NNy;
ndof=5;
xx=aload("wallke.nod.tmp");
x=xx(1:NNy:nnod,1);
if !exist("state") 
  args=[];
else
  args=split(state,"/");
  eval(["clear " state]);
endif

nargs = rows(args);
if nargs<3
  indx=1;
else
  indx = str2num (args(3,:));
endif

if nargs<2
  filename = "save.state";
else 
  filename = deblank(args(2,:));
  if length(filename)==0
    filename="outvector0.out";
  endif
endif

if nargs<1;
  ss="";
else
  ss=deblank(args(1,:));
endif

u=read_state(filename,(nnod+NNy),ndof,indx);
out=NNy:NNy:nnod;
vout=u(out,2);
cntr=(NN-1)*NNy+(1:NNy)';
yc=xx(cntr,2);
uc=u(cntr,2);

U=reshape(u(1:NNy*NN,1),NNy,NN)';
V=reshape(u(1:NNy*NN,2),NNy,NN)';
P=reshape(u(1:NNy*NN,3),NNy,NN)';
K=reshape(u(1:NNy*NN,4),NNy,NN)';
E=reshape(u(1:NNy*NN,5),NNy,NN)';
la=u(NNy*NN+(1:NNy),[1 2]);

min(min(K))
min(min(E))

Q=sum(leftscal(diff(x),xcent(V)));

if exist("ss") && length(ss)>0
  eval([ss ".la = la;"]);
  eval([ss ".U = U;"]);
  eval([ss ".V = V;"]);
  eval([ss ".P = P;"]);
  eval([ss ".K = K;"]);
  eval([ss ".E = E;"]);
  eval([ss ".Q = Q;"]);
  clear ss
endif

if exist(state)
  clear state
endif

u_po=[3.6;
      3.95;
      4.15;
      4.33;
      4.45;
      4.56;
      4.65;
      4.70;
      4.72;
      4.75];
u_po=u_po/max(u_po);
x_po = (1:10)'/10;

return
plot(x,vout);
pause
plot(yc,uc);

#  C_mu=0.09;
#  kappa=0.4;
#  y_wall_plus=50;

#dudy=-(v(2)-v(1))/(x(2)-x(1));
dudy = -cloud(x(1:3)-x(1),1,2)*v(1:3);
nu_eff = nu + C_mu * u(1,5)^2/u(1,6);
tau_w = nu_eff*dudy;
ustar = sqrt(tau_w);
f = v(1)/ustar;
yplus=iwallf(f);
yp = yplus*nu/ustar;
tau_w

## Verif non-linear restriction on k and eps
ustar = v(1)/wallf(y_wall_plus);
ustar^2
kp=ustar^2/sqrt(C_mu);
ep=ustar^4/(von_Karman_cnst*nu*y_wall_plus);

## L es el radio
Q=-2*pi*sum(xcent(v).*diff(x).*xcent(x));
Uav=Q/(pi*L^2);
Re = Uav*2*L/nu

				# Re=4000
u_nk_4000 = [4.15;
	4.65;
	5.55;
	6.22;
	6.65;
	7.;
	7.25;
	7.5;
	7.65;
	7.78;
	7.9;
	8];

u_nk_1_1e5=[4.68;
	    5.2;
	    5.95;
	    6.5;
	    6.9;
	    7.2;
	    7.45;
	    7.6;
	    7.72;
	    7.85;
	    7.92;
	    8];
u_nk_1_1e5 = u_nk_1_1e5/max(u_nk_1_1e5);

u_nk_3_2e6 = [5.48;
	      5.88;
	      6.48;
	      6.93;
	      7.22;
	      7.45;
	      7.6;
	      7.73;
	      7.85;
	      7.9];
u_nk_3_2e6 = u_nk_3_2e6/max(u_nk_3_2e6);

x_nk=[0.15;
      0.32;
      0.79;
      1.6;
      2.4;
      3.2;
      4.;
      4.8;
      5.6;
      6.4;
      7.2;
      8];
x_nk = x_nk / max(x_nk);
plot([0;(1-x+yp)/(1+yp)],[0;v]/max(v),x_nk,u_nk_1_1e5,'o')
