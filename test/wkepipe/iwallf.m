function yplus = iwallf(f)

  f1=5;
  f2 = 2.5*log(30)+5.5;
  yplus  = (f<f1) .* f + (f>=f1 & f<f2) .* exp((f+3.05)/5.0) + ...
      (f>f2) .* exp((f-5.5)/2.5);

endfunction
