#usage: 
function res = fs_eq(X)

  y = fs_fun(X(1));
  res = X(2)-y;

endfunction
