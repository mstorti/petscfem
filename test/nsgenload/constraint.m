##__INSERT_LICENSE__
## $Id: constraint.m,v 1.2 2003/01/08 15:49:04 mstorti Exp $
## This is a simple example showing how the constraints
## are treated inside PETSc-FEM

nedof = 10;
tol = 1e-10;

Q=eye(nedof);

cons = [-1.   2     .1    4       0  9;
	0.1   2     -1    4       0 10;
	0.1   2     -1.1  4       0 10];

for k=1:rows(cons)
  consj = cons(k,:);
  row = zeros(1,nedof);
  for l=1:columns(consj)/2
    coef = consj(1); edof = consj(2); consj(1:2)=[];
    row = row + coef*Q(edof,:);
  endfor
  if (max(abs(row))<tol);  continue; endif
  [bid,edof_elim] = max(abs(row));
  coef = row(edof_elim);
  row(edof_elim) = 0;
  row = row/coef;
  col = Q(:,edof_elim);
  Q(:,edof_elim) = zeros(nedof,1);
  Q = Q - col*row;
endfor

Q
