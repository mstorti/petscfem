%usage: returns nodes on a longitudinal (y,z=cte) column of nodes
function II = indxx (jj,kk)

  global N
  II = (N+1)^2*(kk-1) + (jj-1)*(N+1) + (1:N+1)';

endfunction
