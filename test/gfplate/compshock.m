##usage: U2 = compshock(U1[,du,gamma]);
function U2 = compshock (U1,du,ga);

  if nargin<2; du=0; endif
  if nargin<3; ga=1.4; endif

  ## Compute steady shock wave 
  rho1 = U1(1);
  u1 = U1(2);
  p1 = U1(3);

  ## Values upstream
  c1 = sqrt(ga*p1/rho1);
  M1 = u1/c1;

  ga1 = ga-1;
  C  = (M1+1/(ga*M1))^2/(0.5*M1^2+1/ga1);

  coef = [C/2-1,C/ga1-2/ga,-1/ga^2];
  M2 = sqrt(roots(coef));

  [bid,indx] = min(abs(M2-M1));
  abs(M2(indx)-M1)<1e-10 || error("inconsistent result");

  ## M2 is Mach downstream
  M2(indx) = [];

  ## Values downstream
  u2 = u1*sqrt((0.5+1/M1^2/ga1)/(0.5+1/M2^2/ga1));
  c2 = u2/M2;
  rho2 = u1*rho1/u2;
  p2 = c2^2*rho2/ga;

  U2 = [rho2,u2,p2];

endfunction
