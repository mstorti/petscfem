function [f,df] = wallf(yp)

  c1 = -5*log(5)+5;
  c2=5*log(30)+c1-2.5*log(30);
  f = yp.*(yp<5) + (5.0*log(yp)+c1).*(yp>=5 & yp< 30) + \
      (2.5*log(yp)+c2).*(yp>=30);
  df = (yp<5) + (5.0./yp).*(yp>=5 & yp< 30) + \
      (2.5./yp).*(yp>=30);

#  f = yp.*(yp<30)+ (2.5*log(yp)+5.5).*(yp>=30);

endfunction
