ga = 1.4;
ga1 = ga-1;

M1 = 0.3;
p1 = 1;
rho1 = 1;

c1 = sqrt(ga*p1/rho1);
u1 = M1*c1;

C  = (M1+1/(ga*M1))^2/(0.5*M1^2+1/ga1);

coef = [C/2-1,C/ga1-2/ga,-1/ga^2];
M2 = sqrt(roots(coef));

[bid,indx] = min(abs(M2-M1));
abs(M2(indx)-M1)<1e-10 || error("inconsistent result");

M2(indx) = [];

u2 = u1*sqrt((0.5+1/M1^2/ga1)/(0.5+1/M2^2/ga1));
c2 = u2/M2;
rho2 = u1*rho1/u2;
p2 = c2^2*rho2/ga;

U1 = [rho1,u1,p1]
U2 = [rho2,u2,p2]

if 1
  M = (0.2:0.05:50)';
  
  f = (M.^2+1/ga)./(0.5*M.^2+1/ga1);
  semilogx(M,f);
endif
