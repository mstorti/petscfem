ndof=3;
nstep=1000;
df=0.1;
dt=0.1;


A=rand(ndof);
A=expm(A+A');

C=df*A;

x0=zeros(ndof,1);
f=rand(ndof,1);
II=eye(ndof);

JJ = II/dt^2+C/(2*Dt)+A/2;

xf=x0;
tau = 5/min(eig(A));
beta = 1-exp(-dt/tau);

xo = x0;
xoo= x0;
xh = zeros(nstep+1,ndof);
xfh=xh;
xh(1,:)=xo';
xfh(1,:)=xo';
for k=1:nstep
  x = xo;
  res=(x-2*xo+xoo)/dt^2 + C* (x-xoo)/(2*dt) + A*(x+xoo)/2-f;
  x=x-JJ\res;
  
  xf = beta * x + (1-beta) * xf;

  xoo=xo;
  xo=x;
  xh(k+1,:)=x';
  xfh(k+1,:)=xf';
end
  
