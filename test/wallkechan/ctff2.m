## usage: w = ctff2 (x)
function w = ctff2 (x)

# $Id: ctff2.m,v 1.2 2001/07/02 14:24:19 mstorti Exp $
#__INSERT_LICENSE__

  s=size(x);
  x = vec(x);

  w=zeros(size(x));

#  disp('indx = find(x>=60);');
#    indx = find(x>=60);
#    if length(indx)>0
#      w(indx) = x(indx);
#    endif

#    indx = find(x<=-60);
#  #  disp('indx = find(x<=-60);');;
#    if length(indx)>0
#      w(indx) = zeros(size(indx));
#    endif

#  disp('indx = find(abs(x)<=1e-7);');
  indx = find(abs(x)<=1e-7);
  if length(indx)>0
    xx = x(indx);
    w(indx) = 0.5*exp(xx)./(1+xx.^2/6);
  endif

#  disp('indx = find(abs(x)>1e-7 & abs(x)<60 & x>0);');
  indx = find(abs(x)>1e-7 & x>0);
  if length(indx)>0
    xx = x(indx);
    w(indx) = xx./(1-exp(-2*xx));
  endif

#  disp('indx = find(abs(x)>1e-7 & abs(x)<60 & x<0);');
  indx = find(abs(x)>1e-7 & x<0);
  if length(indx)>0
    xx = x(indx);
    ee = exp(2*xx);
    w(indx) = xx.*ee./(ee-1);
  endif

  w = reshape(w,s(1),s(2));

endfunction
