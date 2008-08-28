source ("data.m.tmp");
  
uref =  [9.9826565720e-01
         9.9979737827e-01
         1.0000435548e+00
         1.0004042038e+00
         1.0004424587e+00
         1.0004574859e+00
         1.0004424587e+00
         1.0004042038e+00
         1.0000435548e+00
         9.9979737827e-01
         9.9826565720e-01];

ur = aload("strip2d.state-right.tmp");
ul = aload("strip2d.state-left.tmp");
erro_r = merr(ur(1:Ny+1,4)-uref);
erro_l = merr(ul(Nx*(Ny+1)+(1:Ny+1),4)-uref);

tol = 1e-7;
printf(["Flow reversal. \n" \
	"Error (flow to right) < tol OK ? %d (error = %g, tol = %g)\n" \
	"Error (flow to left)  < tol OK ? %d (error = %g, tol = %g)\n" \
	"Test OK ? %d\n"], \
	erro_r<tol,erro_r,tol,
	erro_l<tol,erro_l,tol,
	erro_r<tol && erro_l<tol);
