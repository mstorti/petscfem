#usage: 
function f = filana (a,b,N);

  if nargin<3
    N=128;
  endif
  x = zeros(N,1);
  A = x;
  B = x;

  x(N/2)=1;
  A(N/2+[0 1 2]) = a;
  B(N/2+[0 1 2]) = b;

  f = fft(B)./fft(A);

endfunction
