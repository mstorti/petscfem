function res = pp(ustar)

  global visco y_wall u
  yplus = y_wall*ustar/visco;
  res = u - ustar * wallf(yplus);

endfunction 
