%usage: 
function z = myrem (x,y)

#  z  = x - floor(x/y);
  z=rem(x,y);
  z = z + abs(y)*(z<0);

endfunction
