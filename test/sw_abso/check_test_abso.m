source("data.m.tmp");

U=aload("save.state.tmp");

nnodesplag=Nx+2;
h=U(1:Nx,2);
u=U(1:Nx,1);
clear U;

tol_h=4e-3;
tol_u=tol_h;
 
erro_h = sum(abs(diff(h)));
erro_u = sum(abs(diff(u)));

printf("Weak form 0. error_h < tol OK ? %d (error = %g, tol = %g)\n",erro_h<tol_h,erro_h,tol_h);
printf("Weak form 0. error_u < tol OK ? %d (error = %g, tol = %g)\n",erro_u<tol_u,erro_u,tol_u);
