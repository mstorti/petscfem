#usage: 
function x = parabf (y)

  global parab_param
  a = parab_param.a;
  expo = parab_param.expo;
  coe2 = parab_param.coe2;
  x = a*(-1+(y/a).^2/2+coe2*abs(y/a).^expo);

endfunction
