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
  xe = 1.5*rand(1,2)-0.25;			# point to project
  ## xe = [1 1];
  for k=1:nelem1
    nodes = ico1(k,:);
    C = [x1(nodes,:)';
	 ones(1,3)];
    invC = inv(C)';
    invC(ndim,:)=0;
    b = [xe';1];
    flag = zeros(3,1);			# indices for which restrictions are active
    while 1
      indx = find(flag);
      C2 = C;
      if length(indx)>0;
	C2(:,indx) = -invC(:,indx);
      endif
      L = C2\b;
      bad = find(L<0);
      if length(bad)==0; 
	break; 
      else
	flag(bad) = !flag(bad);
      endif
    endwhile
    indx = find(flag);
    L(indx)=0;
    xproj = (C(1:2,:)*L)';
    xprojh=[xprojh;
	    xproj];
    xeh=[xeh;
	 xe];
  endfor
endfor

plot(xeh(:,1),xeh(:,2),'og',xprojh(:,1),xprojh(:,2),'or',\
     [xeh(:,1)';xprojh(:,1)'], [xeh(:,2)';xprojh(:,2)'],b);
