source("~/.octaverc");

u=aload("oscsome4_30deg.sal");
err = norm(u(:,3) - sin((1:50)'/16*2*pi-pi/6));
tol=1e-8;
printf("error = %g,  < %g OK ? %d\n",err,tol,err<tol);
