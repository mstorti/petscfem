#usage: 
function g = fcostg(a0)

  eps = 1e-3;
  f0 = fcost(a0);
  g = zeros(size(a0));
  a = zeros(size(a0));
  for k=1:length(a0)
    a = a0;
    a(k) = a(k)+eps;
    g(k) = (fcost(a)-f0)/eps;
  endfor

endfunction
