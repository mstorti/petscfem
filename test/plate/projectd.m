function xpro = projectd (xini,nor,s,fun);

  tol = 1e-10;
  ## We look for a solution of the form xini + s * nor
  s0=0; s1=s;
  r0 = feval(fun,xini+s0*nor);
  r1 = feval(fun,xini+s1*nor);
  while 1
    drds = (r1-r0)/(s1-s0);
    s2 = s1 - r1/drds;
    x2 = xini+s2*nor;
    r2 = feval(fun,x2);
    if (abs(r2)<tol) ; break ; endif
    s0=s1; r0=r1;
    s1=s2; r1=r2;
  endwhile

  xpro = x2;

endfunction
