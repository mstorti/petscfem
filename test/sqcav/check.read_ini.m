source("data.m.tmp");
u = aload("sqcav.read_ini.state.tmp");

printf("v=p=0 OK ? %d\n",merr(u(:,2:3))==0);
tol = 1e-10;
printf("sum(sum(u)) OK ? %d\n",abs(sum(sum(u(:,1)))-(N+1+u_rini*(N-1)^2))<tol);
