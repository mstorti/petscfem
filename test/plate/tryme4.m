X=[0 0;
   0 1;
   0 2;
   0 3;
   1 0;
   1 1;
   1 2;
   1 3;
   2 0;
   2 1;
   2 2;
   2 3];

x=X(:,1); y=X(:,2);
N=rows(X);

A=[ones(N,1) x x.*x  y y.*x y.*x.*x y.*y y.*y.*x y.*y.*x.*x];
A = leftscal(1./l2(A),A);
nodosA = (1:N)';

AA=A(1,:);
nodosAA=1;
A(1,:)=[];
nodosA(1)=[];

for k=1:N-1
  H = A*AA';
  [bid,indx] = min(max(H'));

  AA=[AA;
      A(indx,:)];
  A(indx,:)=[];
  nodosAA=[nodosAA;
	   nodosA(indx)];
  nodosA(indx)=[];
endfor
