#usage: 
function q = fun (d,p)

  if nargin==1; p=2; endif
  q = (sum(d.^p)*sum(d.^(-p)))^(1/p);

endfunction
