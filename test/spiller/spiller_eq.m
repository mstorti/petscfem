#usage: 
function res = spiller_eq(X)

  y = spiller_fun(X(1));
  res = X(2)-y;

endfunction
