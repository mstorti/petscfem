#usage: 
function ok = mesher_check_corner (x1,x2,tol)
  ok = l2(x1-x2) < tol;
  if !ok
    x1,x2
    error("Corner doesn't match");
  endif
endfunction
