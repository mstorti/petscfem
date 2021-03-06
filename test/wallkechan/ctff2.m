## usage: w = ctff2 (x)
function w = ctff2 (x,tol)

# $Id: ctff2.m,v 1.3 2001/07/03 20:00:40 mstorti Exp $
#__INSERT_LICENSE__

  s=size(x);
  x = vec(x);
  r = x/tol-1;

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
    rr = r(indx);
    w(indx) = 0.5*tol*exp(rr)./(1+rr.^2/6)+tol;
  endif

#  disp('indx = find(abs(x)>1e-7 & abs(x)<60 & x>0);');
  indx = find(x>=1e-7);
  if length(indx)>0
    rr = r(indx);
    w(indx) = (x(indx)-tol)./(1-exp(-2*rr))+tol;
  endif

#  disp('indx = find(abs(x)>1e-7 & abs(x)<60 & x<0);');
  indx = find(x<=-1e-7);
  if length(indx)>0
    rr = r(indx);
    ee = exp(2*rr);
    w(indx) = (x(indx)-tol).*ee./(ee-1)+tol;
  endif

  w = reshape(w,s(1),s(2));

endfunction
