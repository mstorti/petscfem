#usage: 
function y = spillway_fun(x)

  global spillway_data
  s = spillway_data;
  y = s.H1-s.C*x.^s.E;		# Height of spillway w.r.t. flat bottom
  y = choose(x<s.L1,y,s.H2*ones(size(x)));

endfunction
