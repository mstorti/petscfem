##__INSERT_LICENSE__
## $Id: proc.m,v 1.3 2003/01/08 15:49:03 mstorti Exp $
if 0
  petsc_data_name="a_ll_000";
  petscload
  A=getblock(a_ll_000);
else
  petsc_data_name="mat.output";
  petscload
  A=getblock(A);
  AA=getblock(AA);
endif
