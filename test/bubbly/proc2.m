##__INSERT_LICENSE__
## $Id: proc2.m,v 1.2 2003/01/08 15:49:03 mstorti Exp $
petsc_data_name="mat.output";
petscload

indxbb=[2 3 5 6 1 4];
A=getblock(A);
A=A(indxbb,indxbb)
AA=getblock(AA);
AA=AA(indxbb,indxbb)

A-AA
merr(A-AA)
