#usage: 
function res = spillway_eq(X)

  y = spillway_fun(X(1));
  res = X(2)-y;

endfunction
