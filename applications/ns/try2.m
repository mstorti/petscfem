##__INSERT_LICENSE__
## $Id: try2.m,v 1.2 2003/01/08 15:54:25 mstorti Exp $
global visco y_wall u

visco = 0.0000104166666666667;
u= 0.143424000000000;
y_wall=0.01;

ustar = fsolve('pp',sqrt(u*nu/y_wall));
tau = ustar^2;
g = tau/u;
yplus = y_wall*ustar/nu;
[f,fp] = wallf(yplus);
dusdua = 1/(f+yplus*fp);
dga = (2*ustar*dusdua-g)/u;

du =  .01*u;
u = u +du;
ustarn = fsolve('pp',sqrt(nu*u/y_wall));
taun = ustarn^2;
gn = taun/u;

dusdu = (ustarn-ustar)/du;
dg = (gn-g)/du;

