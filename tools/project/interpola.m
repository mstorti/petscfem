if 0
  x1 = aload("static_p_blade.nod");
  ico1 = aload("blade.con");
  
  x2 = aload("patran.nod");
else
  x1 = [0,0,0;
	1,0,0;
	1,1,0;
	0,1,0];
  ico1 = [1 2 3;
	  3 4 1];
  x2 = [0.25 0.25 1;
	0.75 0.75 -1];
endif

nnod1 = rows(x1);
nnod2 = rows(x2);
nelem1 = rows(ico1);

a = x1(ico1(:,2),:)-x(ico1(:,1),:);
b = x1(ico1(:,3),:)-x(ico1(:,1),:);
c = pvec(a,b);
A = [a,b,c];

dx = zeros(nelem1,3);
xreco = zeros(size(dx));
for k=1:nnod2
  for l=1:3
    dx(:,l) =  x2(k,l)-x1(ico1(:,1),l);
  endfor
  L = rowgauss(A,dx,3,1);
  xreco = 0.;
  xreco = xreco + leftscal(L(:,1),a);
  xreco = xreco + leftscal(L(:,2),b);
  xreco = xreco + leftscal(L(:,3),c);
  return
  if 0
    for l=1:2
      indx = find(L(:,l)<0);
      if length(indx)>0
	L(indx,l) = 0;
      endif
    endfor
    L(:,3) = 1-L(:,1)-L(:,2);
    indx = find(L(:,k)<0);
    L(indx,l) = 0;
  endif
endfor
