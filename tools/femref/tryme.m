tol=1e-5;
nindx = find(abs(xnod(:,3)-1)+abs(xnod(:,2))<tol);
h = 2/N;
nindx2 = find(xnod(:,2)>h/2 & xnod(:,2)<h & abs(xnod(:,3)-1)<tol);
nindx = [nindx;nindx2];
nindx = sortby(xnod(nindx,1),nindx);
return

xaxis = xnod(nindx,1);
utheta = xaxis./(dd^2+xaxis.^2);
