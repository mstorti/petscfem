function res = wallf1(ustar)

  global y nu U

  yp = ustar*y/nu;
  yp

  res = wallf(yp)-U/ustar;

endfunction
