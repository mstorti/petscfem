## Copyright (C) 2002 Mario A. Storti
##
## This file is part of Octave.
##__INSERT_LICENSE__
## $Id: check.swturb.m,v 1.1 2002/11/03 01:30:18 mstorti Exp $

##usage: 

## Author: Mario Storti
## Keywords: tests, petscfem, turbulence, channel

uref=aload("swturb.ans");
uref=uref(:,[1 1 1 2 3]);
uref(:,2)=0;
uref(:,3)=0.1;
      
u=aload("swturb.out.tmp");

err_v = max(abs(u-uref));
printf("Max error per component: %f %f %f %f %f\n",err_v);

abs_ref = max(u);
rel_err_v = err_v([1 3 4 5])./abs_ref([1 3 4 5]);
erro = max(err_v);

printf("Max rel. error per component: %f %f %f %f\n",rel_err_v);

tol = 1e-6;
printf("Max rel error OK ? %d [tol %f]\n",max(rel_err_v)<tol,tol);
