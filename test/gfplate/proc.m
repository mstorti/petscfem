## $Id: proc.m,v 1.3 2005/01/09 00:09:25 mstorti Exp $
source("data.m.tmp");

use_previously_loaded_u = 0;
if ! (use_previously_loaded_u && exists("u"))
  u = aload("gfplate.state.tmp");
endif

nnod = (Nx+1)*(Ny+1);
q = -u(nnod+(1:Nx+1),1);
u = u(1:nnod,:);
xnod = aload("gfplate.nod.tmp");
xnod = xnod(1:nnod,:);
nnod = size(xnod,1);
y = xnod(1:Ny+1,2);
x = xnod(1:(Ny+1):nnod,1);

title("density");
xlabel("x");
plot(x,reshape(u(:,1),Ny+1,Nx+1)');
pause;
xlabel("y");
plot(y,reshape(u(:,1),Ny+1,Nx+1));
pause;

title("Ux");
xlabel("x");
plot(x,reshape(u(:,2),Ny+1,Nx+1)');
pause;
xlabel("y");
plot(y,reshape(u(:,2),Ny+1,Nx+1));
pause;

title("Uy");
xlabel("x");
plot(x,reshape(u(:,3),Ny+1,Nx+1)');
pause;
xlabel("y");
plot(y,reshape(u(:,3),Ny+1,Nx+1));
pause;

title("pressure");
xlabel("x");
plot(x,reshape(u(:,4),Ny+1,Nx+1)');
pause;
xlabel("y");
plot(y,reshape(u(:,4),Ny+1,Nx+1));
pause;

title("Temperature");
xlabel("x");
plot(x,reshape(u(:,4)./(u(:,1)*Rgas),Ny+1,Nx+1)');
pause;
xlabel("y");
plot(y,reshape(u(:,4)./(u(:,1)*Rgas),Ny+1,Nx+1));
pause;

indx = find(x<Lplate);
q(indx)=0;

title("Qwall");
xlabel("x");
plot(x,q);
pause;

title("");
xlabel("");
