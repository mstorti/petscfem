function akl = compkl(Ant,k,l);

  ndof=round(sqrt(columns(Ant)));
  nx=round(sqrt(rows(Ant)));
  akl = reshape(Ant(:,(k-1)*ndof+l),nx,nx);

endfunction 
