function flux=upflux(U,n);

  global gravity addvisc

  g=gravity;

  NN=size(U,1);
  Uh=(U(2:NN,:)+U(1:NN-1,:))/2;
  h=Uh(:,1);
  u=Uh(:,2)./h;
  v=Uh(:,3)./h;

  un=u*n(1)+v*n(2);
  flux=[h.*un h.*un.*u+(g*n(1)/2)*h.^2  h.*un.*v+(g*n(2)/2)*h.^2];
  An=[zeros(NN-1,1),-un.*u+g*h*n(1),-un.*v+g*h*n(2),n(1)*ones(NN-1,1),n(1)*u+un,n(1)*v,n(2)*ones(NN-1,1),n(2)*u,n(2)*v+un];

  for k=1:NN-1
    U1=U(k,:);
    U2=U(k+1,:);
    A=reshape(An(k,:),3,3);
    [V,d]=eig(A);
    d=abs(diag(d));
    md=max(d);
    d=d+.3*(md-d);
    aJJ=V*leftscal(d,inv(V));
    flux(k,:) = flux(k,:)  - (aJJ*(U2-U1)')'/2;
  end

