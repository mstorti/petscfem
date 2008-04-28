source ("data.m.tmp");

uref =  [1.1046054554e+00  
         9.6474601775e-01  
         9.9571585206e-01  
         9.9515405613e-01  
         9.9655834530e-01  
         9.9515405613e-01  
         9.9571585206e-01  
         9.6474601775e-01  
         1.1046054554e+00 ];

ur = aload("strip2d.state-right.tmp");
ul = aload("strip2d.state-left.tmp");
return
erro = merr(u(2:Ny,4)-uref);

printf(["Square cavity at Re=1000. " \
	"Error < tol OK ? %d (error = %g, tol = %g)\n"], \
	erro<tol,erro,tol);
