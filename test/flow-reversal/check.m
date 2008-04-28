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
erro_r = merr(ur(2:Ny,4)-uref);
erro_l = merr(ul(Nx*(Ny+1)+(2:Ny),4)-uref);

tol = 1e-7;
printf(["Flow reversal. \n" \
	"Error (flow to right) < tol OK ? %d (error = %g, tol = %g)\n" \
	"Error (flow to left)  < tol OK ? %d (error = %g, tol = %g)\n" \
	"Test OK ? %d\n"], \
	erro_r<tol,erro_r,tol,
	erro_l<tol,erro_l,tol,
	erro_r<tol && erro_l<tol);
