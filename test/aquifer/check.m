##__INSERT_LICENSE__
## $Id: check.m,v 1.2 2003/01/08 15:49:03 mstorti Exp $
source("data.m.tmp");

u=aload(["aqui." casen ".tmp"]);
x=aload("aqui.nod.tmp");

ix=(1:2:2*Nx+2)';
x=x(ix,1);
u=u(ix);
du=diff(u)./diff(x);

xx=xcent(x);
uu=xcent(u);

eta = eta0 + (etaL-eta0)/Lx*xx;

q=du.*(uu-eta);

rele = 2*abs(max(q)-min(q))/abs(max(q)+min(q));
tol=1e-3;

printf("max q: %f, min q: %f, relative deviation: %g, tol: %g\n",
       max(q),min(q),rele,tol);
printf("test OK ? %d\n",rele<tol);
