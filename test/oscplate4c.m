source("~/.octaverc");

u=aload("oscsome4c.sal");
err = norm(u(:,3) - cos((1:50)'/16*2*pi));
tol=1e-8;
printf("error = %g,  < %g OK ? %d\n",err,tol,err<tol);
