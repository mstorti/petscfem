##__INSERT_LICENSE__
## $Id: fcost.m,v 1.1 2003/01/10 12:38:54 mstorti Exp $

function cost = fcost(a,f)
  
  K = [1 1 1 0  0  0 1;
       0 0 0 1  1  1 1;
       0 1 2 0 -1 -2 0;
       0 1 0 0  0  0 a(1);
       0 0 1 0  0  0 a(2);
       0 0 0 0  0  1 a(3)];

  b = K(:,7);
  K(:,7) = [];
  x = K\b;
  a = x(1:3);
  b = x(4:6);
  if nargin==2
    cost = x;
    return;
  endif

  N=128;
  f = filana(a,b,N);
  n1 = round(N/4);
  n2 = round(N/2);
  cost = 50*max(abs(f(n1:n2))) + max(abs(abs(f(1:n1))-1));
  
endfunction 
