source("data.m.tmp");

phi=aload("save.state.tmp");
## phi=aload("advec.ini.tmp");
xx=aload("advec.nod.tmp");
y=xx(1:M+1,2);
x=xx(1:(M+1):rows(xx),1);
phi = reshape(phi,M+1,N+1)';
mesh(phi)
pause
plot(x,phi);
pause
plot(y,phi');
