#usage: 
function [xnod,icone]= mesher_isomap(mesh)

  XNOD = mesh.XNOD;
  ICONE = mesh.ICONE;
  HH = mesh.HH;

  size(HH,2)==1 || size(HH,2)==4 || error("refinement info should have dimension 1 or 4");

  maxnode = sqrt(length(HH));
  edgex = ICONE([1 2]);
  edgex = min(edgex)+maxnode*(max(edgex)-1);
  m = HH(edgex,1);
  x_ref = [1 0 1];
  if size(HH,2)==4; x_ref = HH(edgex,2:4); endif
  edgey = ICONE([1 4]);
  edgey = min(edgey)+maxnode*(max(edgey)-1);
  if x_ref(1)==0; x_ref(1)=1; endif
  if x_ref(3)==0; x_ref(3)=1; endif
  n = HH(edgey,1);
  y_ref = [1 0 1];
  if size(HH,2)==4; y_ref = HH(edgey,2:4); endif
  if y_ref(1)==0; y_ref(1)=1; endif
  if y_ref(3)==0; y_ref(3)=1; endif

  if m==0 || n==0; 
    printf("Not defined refinement!!\n"); 
    ICONE,m,n
    if m==0, m=10; endif
    if n==0, n=10; endif
  endif

  w_not_ref = zhomo([0 1 0 1],m+1,n+1);
  [xinod_nr,icone] = pfcm2fem(w_not_ref);
  clear w_not_ref
  w = zhomo([0 1 0 1],m+1,n+1,[x_ref y_ref]);
  [xinod,icone] = pfcm2fem(w);
  icone = icone(:,[1 4 3 2]);

  xi = onedstr(x_ref,m+1);
  eta = onedstr(y_ref,n+1);
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
  mesher_check_corner(X12(1,  :),X14(1,:),tol);
  X1 = X12(1,  :);
  mesher_check_corner(X12(m+1,:),X23(1,:),tol);
  X2 = X12(m+1,  :);
  mesher_check_corner(X23(n+1,:),X43(m+1,:),tol);
  X3 = X23(n+1,  :);
  mesher_check_corner(X14(n+1,:),X43(1,:),tol);
  X4 = X14(n+1,  :);

  xnod = zeros(size(xinod));
  xi  = xinod(:,1);
  eta = xinod(:,2);
  ii = 1+round(xinod_nr(:,1)*m);
  jj = 1+round(xinod_nr(:,2)*n);
  XBL = \
      (1-xi).*(1-eta)*X1+ \
      xi.*(1-eta)*X2+     \
      xi.*eta*X3 +        \
      (1-xi).*eta*X4;
  XS = (1-xi)*X1+xi*X2;
  XN = (1-xi)*X4+xi*X3;
  XW = (1-eta)*X1+eta*X4;
  XE = (1-eta)*X2+eta*X3;
  XPL = \
      leftscal((1-eta),(X12(ii,:)-XS)) + leftscal(eta,(X43(ii,:)-XN)) + \
      leftscal((1-xi),(X14(jj,:)-XW)) + leftscal(xi,(X23(jj,:)-XE));
  xnod = XBL + XPL;

endfunction
