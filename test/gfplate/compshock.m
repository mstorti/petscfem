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

  ## Solution is found by considering that 
  ## (Ma+1/(ga*Ma))^2/(0.5*Ma^2+1/ga1) is constant
  ## across the shock (Rankine-Hugoniot relations)
  ga1 = ga-1;
  C  = (M1+1/(ga*M1))^2/(0.5*M1^2+1/ga1);

  ## Solves polynomial equation 
  coef = [C/2-1,C/ga1-2/ga,-1/ga^2];
  M2 = sqrt(roots(coef));

  ## This should have two solutions, one
  ## is M1 and the other should be M2
  [bid,indx] = min(abs(M2-M1));
  abs(M2(indx)-M1)<1e-10 || error("inconsistent result");

  ## M2 is Mach downstream
  M2(indx) = [];

  ## Take care that, for a given `gamma' there is a minimum
  ## downstream Mach possible (1/sqrt(7)=0.37796 for gamma=1.4).
  ## If the 'M1' is lower than this value, then you can get complex
  ## solutions to the equation. 
  abs(imag(M2))<1e-10 || error("can't find a real solution");

  ## Values downstream
  u2 = u1*sqrt((0.5+1/M1^2/ga1)/(0.5+1/M2^2/ga1));
  c2 = u2/M2;
  rho2 = u1*rho1/u2;
  p2 = c2^2*rho2/ga;

  ## State downstream
  U2 = [rho2,u2,p2];

endfunction
