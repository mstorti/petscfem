### key proc4.m
### $Id: proc4.m,v 1.1 2005/01/21 17:10:47 mstorti Exp $

source("data.m.tmp");
U=aload("gfabso.state.tmp");

rho = reshape(U(:,1),Nx+1,2);
rho=rho(1:Nx-1,:);

u = reshape(U(:,2),Nx+1,2);
u = u(1:Nx-1,:);

ndof=4;
uu=U([1:Nx-1,Nx+1],:);
gasdata.gamma = gamma;
uri = primi2ri(uu,gasdata);
