##__INSERT_LICENSE__
## $Id: trapez.m,v 1.2 2003/01/08 15:49:03 mstorti Exp $
nelem=10;
h=1./nelem;
A=(1/h)*(2*diag(ones(nelem-1,1))-diag(ones(nelem-2,1),-1)-diag(ones(nelem-2,1),+1));

alpha=.5;
Dt=0.5;
N=round(5/Dt);
f0=zeros(nelem-1,1); f0(1)=1/h;
x0=zeros(nelem-1,1);
omega=pi/2;

K = (eye(size(A))/Dt+alpha*A);
iK=inv(K);

a=(i*omega*eye(size(A))+A)\f0;
b=x0-a;
[V,D]=eig(A);
iV=inv(V);
D=diag(D);

x=x0;
t=0;
xh=[x'];
xah=[x0'];
for j=0:N
  ta = t+alpha*Dt;
  dx = iK*(f0*cos(omega*ta)-A*x);
  x=x+dx;
  xh=[xh ; x'];

  t = t+Dt;
  % analitico
  xho=a*exp(i*omega*t);
  xnho=V*diag(exp(-D*t))*iV*b;

  xah=[xah; real(xho+xnho)'];
endfor

error=sqrt(sum((xah-xh).^2)*Dt)
