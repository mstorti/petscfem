## usage: 
##    smoothes the solution for the cubic cavity problem
##    on a structured mesh
function us = pf_smooth (u,smooth,nstep)

  if nargin<3, nstep=1;   endif 
  if nargin<2, smooth=1; endif
    
  n = rows(u)^(1/3)-1;
  tol = 1e-10;
  abs(round(n)-n) < tol || error("Not a structured cubic mesh");
  n = round(n);

  ## Nodes to the Nord, South, East, West, Up, Down
  N = (n+1)^3;
  n1 = n+1;

  here = (0:N-1)';
  JJ = floor(here/n1);
  KK = floor(JJ/n1);
  JJ = rem(JJ,n1);
  II = rem(here,n1);

  us = u;

  coef = smooth/12;
  for k=1:nstep
    uss = (1-smooth/2)*us;

    neighbor = choose(II!=n,here+1,here-1);
    uss = uss + coef*us(neighbor+1,:);
    neighbor = choose(II!=0,here-1,here+1);
    uss = uss + coef*us(neighbor+1,:);

    neighbor = choose(JJ!=n,here+n1,here-n1);
    uss = uss + coef*us(neighbor+1,:);
    neighbor = choose(JJ!=0,here-n1,here+n1);
    uss = uss + coef*us(neighbor+1,:);

    neighbor = choose(KK!=n,here+n1^2,here-n1^2);
    uss = uss + coef*us(neighbor+1,:);
    neighbor = choose(KK!=0,here-n1^2,here+n1^2);
    uss = uss + coef*us(neighbor+1,:);
    
    us = uss;
  endfor

endfunction
