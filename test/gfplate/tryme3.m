## $Id: tryme3.m,v 1.1 2005/02/19 15:45:14 mstorti Exp $
ga = 1.4;
ga1 = ga-1;

M1 = 0.4;
rho1=1;
p1=1/ga;

c1 = sqrt(ga*p1/rho1);
u1 = c1*M1;
U1 = [rho1,u1,p1];

M2 = (0.1:0.01:20)';

## Get u from energy conservation and mach
ht_1 = 0.5*u1^2+c1^2/ga1;
u2 = sqrt(ht_1./(0.5+1./(M2.^2*ga1)));
c2 = u2./M2;
rho2 = rho1*u1./u2;
p2 = c2.^2.*rho2/ga;


flux1 = [rho1*u1, \
	 rho1*u1^2+p1, \
	 (0.5*u1^2+c1^2/ga1)];

flux2 = [rho2.*u2, \
	 rho2.*u2.^2+p2, \
	 (0.5*u2.^2+c2.^2/ga1)];

merr(flux2(:,1)-flux1(:,1)) < 1e-10 || \
    error("mass balance not satisfied.");

merr(flux2(:,3)-flux1(:,3)) < 1e-10 || \
    error("energy balance not satisfied.");

plot(M2,flux2(:,2)-flux1(:,2));
