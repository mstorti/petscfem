%usage: returns nodes on a vertical (x,z=cte) column of nodes
function JJ = indxy (ii,kk)

  global N
  JJ = (N+1)^2*(kk-1) + (0:N)'*(N+1) + ii;

endfunction
