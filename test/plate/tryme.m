m0=1; xi0=0. ; w0=1;
kc = w0^2*0.2;

A0=[0 1;
    -w0^2 -xi0*w0];

A=zeros(4);

A(1:2,1:2)=A0;
A(3:4,3:4)=A0;

ndof = 4;
Ac=zeros(ndof);
Ac(2,3) = kc;
Ac(4,1) = kc;

x=[0 1 0 0]';
Dt=2*pi/w0/32;

Gm = (eye(ndof)/Dt-A/2);
Gp = (eye(ndof)/Dt+A/2);

Xh = x';
xold = x;
nold = 5;
xold = xold(:,ones(nold,1));

nstep = 100;
for k=1:nstep
  xstar = 2*xold(:,1)-xold(:,2);
  x = Gm \ (Gp*x+Ac*(xstar+x));
  Xh = [Xh;
	x'];
  xold = [x xold(:,1:nold-1)];
endfor
