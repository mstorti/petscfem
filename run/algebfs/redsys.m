function y = redsys (x);

  global AA A B C D H bH
  n=rows(D);
  if rows(x)!= n-1; 
    disp("error: non conforming vector"); 
  endif
  
  x=[x;
     -sum(x(1:n-1))];
  
  y = H*x;
  y = [y(1:n-1)]-bH;

endfunction 
