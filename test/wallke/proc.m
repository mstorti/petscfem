source("data.m.tmp");
NN=2*(N+1);
x=aload("wallke.nod.tmp");
x=x(1:2:NN,1);
u=aload("save.state");
v=u(1:2:NN,2);
plot(x,v)

Chi=0.4;
C_mu=0.09;

dudy=(v(2)-v(1))/(x(2)-x(1));
nu_eff = nu + C_mu * u(1,4)^2/u(1,5);
tau_w = nu_eff*dudy;
ustar = sqrt(tau_w);
yplus = v(1)/ustar;
yp=nu*yplus/ustar;
kp=ustar^2/sqrt(0.09);
ep=ustar^3/Chi/yp;

Q=sum(xcent(v).*diff(x))
Uav=Q/L
Re = Uav*2*L/nu
