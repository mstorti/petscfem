function piecewtanh(slope,file,n);
  if nargin<3; n=20; endif
  x=(0:n)'/n;
  y=tanh(slope*x)/tanh(slope);
  asave(file,[-10 0;x y; 10 1]);
endfunction 
