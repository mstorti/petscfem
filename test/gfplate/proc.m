## $Id: proc.m,v 1.1 2005/01/08 14:34:26 mstorti Exp $
source("data.m.tmp");

u=aload("gfplate.state.tmp");
xnod = aload("gfplate.nod.tmp");
nnod = size(xnod,1);
y = xnod(1:Ny+1,2);
x = xnod(1:(Ny+1):nnod,1);

title("density");
plot(x,reshape(u(:,1),Ny+1,Nx+1)');
pause;
plot(y,reshape(u(:,1),Ny+1,Nx+1));
pause;

title("Ux");
plot(x,reshape(u(:,2),Ny+1,Nx+1)');
pause;
plot(y,reshape(u(:,2),Ny+1,Nx+1));
pause;

title("Uy");
plot(x,reshape(u(:,3),Ny+1,Nx+1)');
pause;
plot(y,reshape(u(:,3),Ny+1,Nx+1));
pause;

title("pressure");
plot(x,reshape(u(:,4),Ny+1,Nx+1)');
pause;
plot(y,reshape(u(:,4),Ny+1,Nx+1));
pause;

title("Temperature");
plot(x,reshape(u(:,4)./(u(:,1)*Rgas),Ny+1,Nx+1)');
pause;
plot(y,reshape(u(:,4)./(u(:,1)*Rgas),Ny+1,Nx+1));
pause;

title("");
