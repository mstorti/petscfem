## Copyright (C) 2002 Mario A. Storti
##
## This file is part of Octave.
##__INSERT_LICENSE__
## $Id: checknutcyl.m,v 1.2 2002/08/19 16:06:15 mstorti Exp $

## Author: Mario Storti
## Keywords: petscfem-test, nutating-cylinder

f0=aload("cylinder.force.nuta0.tmp");
f15=aload("cylinder.force.nuta15.tmp");

Mz = -(f15(:,6)-f0(:,6))*950*0.738; # In lbf ft
Mz_ref = 1.217;			# PETSc-FEM returns this for this
				# mesh (the tru value is around 1.55). 
tol = 0.05;			# relative

rel_err = abs(Mz/Mz_ref-1);
printf("Roll moment OK ? %d, Mz = %f\n",rel_err<tol,Mz);
printf("Mz_ref %f, rel.error: %f, rel.error.tol. %f\n",Mz_ref,rel_err,tol);
       
