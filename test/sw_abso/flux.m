#usage: 
function F = flux (u,h)

  global g
  F = [u^2*h+g*h^2/2; u*h];
  Fstar = [(u/2)*(5/2*g*h+0.5*u^2)*sqrt(h/g); \
	   (0.5*u^2+g*h)*(2/3)*sqrt(h/g)];
  F = (F+Fstar)/2;

endfunction
