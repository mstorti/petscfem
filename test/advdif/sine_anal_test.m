sine

x=[1 1];
source("sine.data");
source("sine.data.res");

uu=aload("sine.some");

ky=pi/2/Ly;

a=D;
b=-u;
c=-(D*ky^2+i*omega);

s=sqrt(b^2-4*a*c);
la1 = (-b + s)/2/a;
la2 = (-b - s)/2/a;

A=[1 1;
   la1*exp(la1*Ly) la2*exp(la2*Ly)];
ab = A\[1;0];

Nx=nx+1;
Ny=ny+1;

N=Nx*Ny;
Ix=(1:ny+1:N)';
xx=xnod(Ix,1);

Iy=(1:Ny+1)';
yy=xnod(Iy,2);

uuh=aload("save.state");
uc=ab(1)*exp(la1*xx)+ab(2)*exp(la2*xx);
phase=exp(i*omega*nstep*Dt/T);
#plot(xx,[real(uc*phase) uuh(Ix)]);
#plot(xx,[imag(uc*phase) uuh(ny+Ix)]);
#pause

uh=aload("sine.some");
uh=uh(:,2);
uuc=sin(ky*Ly)*(ab(1)*exp(la1*x(1))+ab(2)*exp(la2*x(1)));
t=(1:nstep)'*Dt;
uana=imag(uuc*exp(i*omega*t));
#plot(t,[uh uana])

du=uh-uana;
## Error over last period
n1=nstep-N_step_period+1;

fid=fopen("sine_anal_test.sal","a");
fprintf(fid,"Dt= %e   error= %e\n",Dt,norm(du(n1:nstep))*Dt);
fclose(fid);

