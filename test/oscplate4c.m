##__INSERT_LICENSE__
## $Id: oscplate4c.m,v 1.2 2003/01/08 15:49:03 mstorti Exp $
source("~/.octaverc");

u=aload("oscsome4c.sal");
err = norm(u(:,3) - cos((1:50)'/16*2*pi));
tol=1e-8;
printf("error = %g,  < %g OK ? %d\n",err,tol,err<tol);
