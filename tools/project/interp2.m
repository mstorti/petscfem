if 0
  x1 = aload("static_p_blade.nod");
  ico1 = aload("blade.con");
  
  x2 = aload("patran.nod");
else
  x1 = [0,0;
	1,0;
	1,1;
	0,1];
#   ico1 = [2 4 1;
# 	  4 2 3];
  ico1 = [2 4 1];
  x2 = [0.25 0.25;
	0.75 0.75];
endif

ndim = 3;
nnod1 = rows(x1);
nnod2 = rows(x2);
nelem1 = rows(ico1);

C = zeros(ndim);
C2 = C;
xeh=[];
xprojh=[];
for j=1:100
  xe = 3*rand(1,2)-1;			# point to project
  for k=1:nelem1
    nodes = ico1(k,:);
    C = [x1(nodes,:)';
	 ones(1,3)];
    invC = inv(C)';
    invC(ndim,:)=0;
    b = [xe';1];
    L = C\b;
    indx = find(L<0);
    if length(indx) < 0;
      xproj = xe;
    else
      C2 = C;
      C(:,indx) = invC(:,indx);
      L2 = C2\b;
      L2(indx) = 0;
      xproj = (C(1:2,:)*L2)';
    endif
    xprojh=[xprojh;
	    xproj];
    xeh=[xeh;
	 xe];
  endfor
endfor

plot(xeh(:,1),xeh(:,2),'og',xprojh(:,1),xprojh(:,2),'or',\
     [xeh(:,1)';xprojh(:,1)'], [xeh(:,2)';xprojh(:,2)'],b);
