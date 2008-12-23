###key checkinv.m
### $Id: $

source("data.m.tmp");

u0 = aload("gfabso.state-ALE0.tmp");
u0 = u0(1:Nx+1,:);
u1 = aload("gfabso.state-ALE1.tmp");
u1 = u1(1:Nx+1,:);

u1(:,2) = u1(:,2) + uref;

tol = 0.01;
erro = merr(u1-u0);
printf("error < tol OK? %d, (error %f, tol %f)\n",erro<tol,erro,tol);
