function [f,df] = wallf(yp)

  f = yp.*(yp<5) + (5.0*log(yp)-3.05).*(yp>=5 & yp< 30) + \
  (2.5*log(yp)+5.5).*(yp>=30);
  df = (yp<5) + (5.0./yp).*(yp>=5 & yp< 30) + \
      (2.5./yp).*(yp>=30);

#  f = yp.*(yp<30)+ (2.5*log(yp)+5.5).*(yp>=30);

endfunction
