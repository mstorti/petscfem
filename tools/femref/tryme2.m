## $Id: tryme2.m,v 1.1 2005/01/13 20:25:03 mstorti Exp $

surf_con = aload("cube.surf-con-red.tmp")+1;
xs = aload("cube.surf-nod.tmp");
gusurf = aload("cube.grad-un.tmp");
wsurf = [+gusurf(:,10)-gusurf(:,7), \
	 -gusurf(:,3)+gusurf(:,9), \
	 +gusurf(:,5)-gusurf(:,2)];
awsurf = l2(wsurf);

tol=1e-5;
indx = find(xs(:,3)>1-tol);

rho = l2(xs(:,1:2));

NN=100;
x = (0:NN)'/NN*max(rho);
wa = 2*dd^2./(dd^2+x.^2).^2;

plot(rho(indx),awsurf(indx),'o',x,wa);
