function Uproy= absorb (Ustar,Uprev,n)

  global gravity

  h=Uprev(1);
  u=[Uprev(2) Uprev(3)]/h;
  un=sum(u.*n);

  A=[0 n;
     -un*u'+gravity*h*n' u'*n+un*eye(2)];
  [v,d]=eig(A);
  Pi=v*diag(1-sign(diag(d)))/2*inv(v);
  Uproy=Uprev+(Pi*(Ustar-Uprev)')';
