## $Id: interp3d.m,v 1.2 2005/02/24 22:30:17 mstorti Exp $

if 0
  xnod1 = aload("static_p_blade.nod");
  ico1 = aload("blade.con");
  
  xnod2 = aload("patran.nod");
else
  xnod1 = [0,0,2;
	   1,0,2;
	   1,1,2;
	   0,1,2];
  ico1 = [2 4 1];
  xnod2 = [0.25 0.25 1;
	   0.75 0.75 1];
endif

ndim = 3;
nnod1 = rows(xnod1);
nnod2 = rows(xnod2);
nelem1 = rows(ico1);

## Normals to panels
a = xnod1(ico1(:,2),:)-xnod1(ico1(:,1),:);
b = xnod1(ico1(:,3),:)-xnod1(ico1(:,1),:);
nor = pvec(a,b);
clear a b

C = zeros(ndim+1);
C2 = C;
xeh=[];
xprojh=[];
maxit=0;
tol = 1e-6;
for j=1:1
  ##  xe = 1.5*rand(1,3)-0.25;			# point to project
  xe = [0.25,1.4,1];
  for k=1:nelem1
    nodes = ico1(k,:);
    C = [xnod1(nodes,:)';
	 ones(1,3)];
    C = [C,[nor(k,:)';0]];
    invC = inv(C)';
    invC(ndim+1,:)=0;
    b = [xe';1];
    flag = zeros(3,1);			# indices for which restrictions are active
    iters=0;
    while 1
      iters = iters+1;
      indx = find(flag);
      C2 = C;
      if length(indx)>0;
	C2(:,indx) = -invC(:,indx);
      endif
      L = C2\b;
      bad = find(L(1:ndim)<-tol);
      if length(bad)==0; 
	break; 
      else
	flag(bad) = !flag(bad);
      endif
      L(indx)=0;
      L(ndim+1)=0;
      xproj = (C(1:ndim,:)*LL)';
      xproj
    endwhile
    maxit = max([maxit iters]);
    indx = find(flag);
    L(indx)=0;
    L(ndim+1)=0;
    xproj = (C(1:ndim,:)*L)';
    xprojh=[xprojh;
	    xproj];
    xeh=[xeh;
	 xe];
  endfor
endfor

printf("maximum %d iters\n",maxit);

plot(xeh(:,1),xeh(:,2),'og',xprojh(:,1),xprojh(:,2),'or',\
     [xeh(:,1)';xprojh(:,1)'], [xeh(:,2)';xprojh(:,2)'],b);
