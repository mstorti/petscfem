function xc=xcent(x)

  n=size(x,1);
  xc=(x(1:n-1,:)+x(2:n,:))/2;