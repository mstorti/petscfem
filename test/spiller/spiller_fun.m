#usage: 
function y = spiller_fun(x)

  global spiller_data
  s = spiller_data;
  y = s.H1-s.C*x.^s.E;		# Height of spiller w.r.t. flat bottom
  y = choose(x<s.L1,y,s.H2*ones(size(x)));

endfunction
