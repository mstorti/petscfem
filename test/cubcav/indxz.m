%usage: returns nodes on a longitudinal (xmy=cte) column of nodes
function KK = indxz (ii,jj)

  global N
  KK = (N+1)^2*(0:N)' + (jj-1)*(N+1) + ii;

endfunction
