#usage: 
function [xnod,icone]= isomap(XNOD,ICONE,m,n)

  w = zhomo([0 1 0 1],m+1,n+1);
  [xinod,icone] = pfcm2fem(w);

  xi = (0:m)'/m;
  eta = (0:n)'/n;
  X12 = zeros(m+1,2);
  X43 = zeros(m+1,2);
  for k=1:m+1
    X12(k,:) = mapbou(ICONE(1),ICONE(2),xi(k),XNOD);
    X43(k,:) = mapbou(ICONE(4),ICONE(3),xi(k),XNOD);
  endfor
  
  X14 = zeros(n+1,2);
  X23 = zeros(n+1,2);
  for k=1:n+1
    X14(k,:) = mapbou(ICONE(1),ICONE(4),eta(k),XNOD);
    X23(k,:) = mapbou(ICONE(2),ICONE(3),eta(k),XNOD);
  endfor
  
  xnod = zeros(size(xinod));
  for k=1:rows(xnod);
    xi  = xinod(k,1);
    eta = xinod(k,2);
    ii = 1+round(xi*m);
    jj = 1+round(eta*n);
    xnod = (

endfunction
