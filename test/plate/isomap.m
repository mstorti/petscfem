#usage: 
function [xnod,icone]= isomap(XNOD,ICONE,HH)

  maxnode = sqrt(length(HH));
  edgex = ICONE([1 2]);
  edgex = min(edgex)+maxnode*(max(edgex)-1);
  m = HH(edgex);
  edgey = ICONE([1 4]);
  edgey = min(edgey)+maxnode*(max(edgey)-1);
  n = HH(edgey);

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
  
  tol = 1e-10;
  merr(X12(1,  :) - X14(1,:))   < tol || error("Corner doesn't match");
  X1 = X12(1,  :);
  merr(X12(m+1,:) - X23(1,:))   < tol || error("Corner doesn't match");
  X2 = X12(m+1,  :);
  merr(X23(n+1,:) - X43(m+1,:)) < tol || error("Corner doesn't match");
  X3 = X23(n+1,  :);
  merr(X14(n+1,:) - X43(1,:))   < tol || error("Corner doesn't match");
  X4 = X14(n+1,  :);

  xnod = zeros(size(xinod));
  for k=1:rows(xnod);
    xi  = xinod(k,1);
    eta = xinod(k,2);
    ii = 1+round(xi*m);
    jj = 1+round(eta*n);
    XBL = \
	X1*(1-xi)*(1-eta)+ \
	X2*    xi*(1-eta)+ \
	X3*    xi*   eta + \
	X4*(1-xi)*   eta;
    XS = X1*(1-xi)+X2*xi;
    XN = X4*(1-xi)+X3*xi;
    XW = X1*(1-eta)+X4*eta;
    XE = X2*(1-eta)+X3*eta;
    XPL = \
	(1-eta)*(X12(ii,:)-XS) + eta*(X43(ii,:)-XN) + \
	(1-xi)*(X14(jj,:)-XW) + xi*(X23(jj,:)-XE);
    xnod(k,:) = XBL + XPL;
  endfor

endfunction
