##__INSERT_LICENSE__
## $Id: proc.m,v 1.7 2003/01/08 15:49:05 mstorti Exp $
source("data.m.tmp");
NN=2*(N+1);
x=aload("wallke.nod.tmp");
x=x(1:2:NN,1);
u=aload("save.state");
v=u(1:2:NN,2);
plot(x,v)

Chi=0.4;
C_mu=0.09;

#dudy=(v(2)-v(1))/(x(2)-x(1));	# first order
dudy = cloud(x(1:3)-x(1),1,2)*v(1:3);  # 2nd order precision
nu_eff = nu + C_mu * u(1,4)^2/u(1,5);
tau_w = nu_eff*dudy;
ustar = sqrt(tau_w);
f = v(1)/ustar;
yplus = iwallf(f);
yp=yplus*nu/ustar;
kp=ustar^2/sqrt(0.09);
ep=ustar^3/Chi/yp;

Q=2*sum(xcent(v).*diff(x))
Uav=Q/(2*L)
Re = Uav*2*L/nu
#plot([0;x;2-x(NN/2:-1:1);2],[0;v;v(NN/2:-1:1);0]/max(v))
#plot([0;x],[0;v]/max(v))

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

u_co=[2.6;
      2.85;
      3.;
      3.12;
      3.3;
      3.43;
      3.52;
      3.65;
      3.85;
      3.97]/3.97;
x_co=[.2;
      .45;
      .63;
      .89;
      1.25;
      1.58;
      1.95;
      2.61;
      3.3;
      4.02]/4.02;

if strcmp(ocase,"poiseuille")
  plot([0;x],[[0;v/max(v)], [0;v]/max(v)],x_po,u_po,"o")
elseif strcmp(ocase,"couette")
  plot([0;x],[[0;v/max(v)], [0;v]/max(v)],x_co,u_co,"o")
endif

