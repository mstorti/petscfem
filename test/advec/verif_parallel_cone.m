##__INSERT_LICENSE__
## $Id: verif_parallel_cone.m,v 1.3 2003/01/08 15:49:03 mstorti Exp $
###key verif.parallel_cone.m

proc

## Error in x position
err_x = abs((xc(1)-xc0(1))-(ux*Tend))/traversed_length;
err_y = abs((xc(2)-xc0(2))-(uy*Tend))/traversed_length;

tol = 4e-2;
printf("Error in x position < tol OK ? %d (err_x %f, tol %f)\n",
       err_x<tol, err_x, tol);
tol = 1e-4;
printf("Error in y position < tol OK ? %d (err_y %f, tol %f)\n",
       err_y<tol, err_y, tol);

max_rate = 4e-3;
xrate = (sx/sx0-1)/traversed_length;
yrate = (sy/sy0-1)/traversed_length;
printf("Max grouth in x second moment OK ? %d (rate %f, max_rate %f)\n",
       xrate<max_rate,xrate,max_rate);
printf("Max grouth in y second moment OK ? %d (rate %f, max_rate %f)\n",
       yrate<max_rate,yrate,max_rate);

