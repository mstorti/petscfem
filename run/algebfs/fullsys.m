function y = fullsys (x);

  global AA A B C D b
  n=rows(D);
  if rows(x)!= 3*n-1; 
    disp("error: non conforming vector"); 
  endif
  
  x=[x;
     -sum(x(2*n+(1:n-1)))];
  
  y = AA*x;
  y = [y(1:3*n-1)]-b;

endfunction 
