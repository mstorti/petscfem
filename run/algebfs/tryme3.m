x=0*b;
omega=0.9;
da=diag(diag(A));

x = gmres(A,diag(A),5,1e-10);

#  rh=[];
#  while 1
#    r=A*x-b;
#    disp(sprintf('%e',norm(r,1)))
#    rh=[rh;norm(r,1)];
#    x=x-omega*(da\r);
#  end
