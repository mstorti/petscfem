##__INSERT_LICENSE__
## $Id: proc.m,v 1.3 2003/01/08 15:49:05 mstorti Exp $
source("data.m.tmp");
NN=2*(N+1);
x=aload("wallke.nod.tmp");
x=x(1:2:NN,1);
u=aload("save.state");
v=u(1:2:NN,2);

#  C_mu=0.09;
#  kappa=0.4;
#  y_wall_plus=50;

dudy=-(v(2)-v(1))/(x(2)-x(1));
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

plot([0;1-x],[0;v]/max(v))

