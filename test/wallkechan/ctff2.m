## usage: w = ctff2 (x)
function w = ctff2 (x)

# $Id: ctff2.m,v 1.1 2001/07/01 23:10:13 mstorti Exp $
#__INSERT_LICENSE__

  s=size(x);
  x = vec(x);

  w=zeros(size(x));

  indx = find(x>=60);
  keyboard
  if length(indx)>0
    w(indx) = x(indx);
  endif

  indx = find(x<=-60);
  keyboard
  if length(indx)>0
    w(indx) = zeros(size(indx));
  endif

  indx = find(abs(x)<=1e-7);
  keyboard
  if length(indx)>0
    xx = x(indx);
    w(indx) = exp(xx)./(1+xx.^2/6);
  endif

  indx = find(abs(x)>1e-7 && abs(x)<60 && x>0);
  keyboard
  if length(indx)>0
    xx = x(indx);
    w(indx) = xx./(1-exp(-2*xx));
  endif

  indx = find(abs(x)>1e-7 && abs(x)<60 && x<0);
  keyboard
  if length(indx)>0
    xx = x(indx);
    ee = exp(2*xx);
    w(indx) = xx.*ee./(ee-1);
  endif

  w = reshape(w,s(1),s(2));

endfunction
