##__INSERT_LICENSE__
## $Id: proc.m,v 1.1 2003/01/08 15:49:04 mstorti Exp $

source("data.m.tmp");

nod_axis = aload("nozzle.fixa_axis.tmp");
nod_axis = [1;
	    nod_axis(:,1)];

xnod = aload("nozzle.nod.tmp");
x = xnod(nod_axis,1);

U = aload("nozzle.state.tmp");
U = U(nod_axis,:);
#plot(x,U(nod_axis,:));

rho = U(:,1);
p = U(:,4);
u2 = U(:,2).^2+U(:,3).^2;
s = p./rho.^ga;
e = p./rho/(ga-1);
h = p./rho*ga/(ga-1)+0.5*u2;
hi = p./rho*ga/(ga-1);
c = sqrt(ga*(ga-1)*e);
M = sqrt(u2)./c;
