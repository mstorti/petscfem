source("~/.octaverc");

tol=1e-10;
u=aload("state_save.sal");
indx=[61:120 61]';

err = norm(u(indx,1) - u(indx(61:-1:1),1));
printf("error in u symmetry = %g,  < %g OK ? %d\n",err,tol,err<tol);

err = norm(u(indx,2) + u(indx(61:-1:1),2));
printf("error in v symmetry = %g,  < %g OK ? %d\n",err,tol,err<tol);

err = norm(u(indx,3) - u(indx(61:-1:1),3));
printf("error in h symmetry = %g,  < %g OK ? %d\n",err,tol,err<tol);

