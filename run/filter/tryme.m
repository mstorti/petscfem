ndof=3;
nstep=1000;
df=0.01;
cc=2;
co=.3;


A=cc*(2*rand(ndof)-1);
A=expm(A+A');
omega = sqrt(max(eig(A)));
dt = co*2*pi/omega;

C=df*A;

x0=zeros(ndof,1);
ff=rand(ndof,1);
II=eye(ndof);
xex=A\ff;

JJ = II/dt^2+C/(2*dt)+A/2;

xf=x0;
xf2=x0;
xfsq=x0;
xfcb=x0;
xav=x0;
xav2=x0;
xav3=x0;
nav=0;
tau = 2/min(eig(A));
beta = 1-exp(-dt/tau);
beta2 = 1-exp(-2*dt/tau);
betaf=0.2;

f=x0;

xo = x0;
xoo= x0;
nmeth=8;
xh = zeros(nstep+1,ndof*nmeth);
xh(1,:)=kron(ones(1,nmeth),xo');
for k=1:nstep

  f = beta * ff + (1-beta) *f;

  x = xo;
  res=(x-2*xo+xoo)/dt^2 + C* (x-xoo)/(2*dt) + A*(x+xoo)/2-f;
  x=x-JJ\res;
  
  xf = beta * x + (1-beta) * xf;
  xfsq = beta * xf + (1-beta) * xfsq;
  xfcb = beta * xfsq + (1-beta) * xfcb;
  xf2 = beta2 * x + (1-beta2) * xf2;
  xav = (nav*xav+x)/(nav+1);
  xav2 = (nav*xav2+xav)/(nav+1);
  xav3 = (nav*xav3+xav2)/(nav+1);
  nav = nav+1;
  xoo=xo;
  xo=x;
  xh(k+1,:)=[x' xf' xf2' xfsq' xfcb' xav' xav2' xav3'];
end

eh = zeros(nstep+1,nmeth);
for lm=1:nmeth
  for ldof=1:ndof
    indx=ndof*(lm-1)+ldof;
    eh(:,lm) = eh(:,lm)+(xh(:,indx)-xex(ldof)).^2;
  end
end
eh=sqrt(eh);
  
