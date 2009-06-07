###key proc.m
### $Id: $

source("mat.output");

function Acond = condense(A);

  if 0
    i1 = [1,3];
    i2 = [2,4];
    A11 = A(i1,i1);
    A12 = A(i1,i2);
    A21 = A(i2,i1);
    A22 = A(i2,i2);
    tol = 1e-6;
    merr(A11-A22)<tol || error("A11!=A22");
    merr(A12-A21)<tol || error("A12!=A21");
    merr(A12-0.5*A11)<tol || error("A12!=0.5*A11");
    Acond = A11+A12+A21+A22;
  else
    Acond = A;
  endif
  
endfunction 

Aana = condense(full(Mat_1));
Afdj = condense(full(Mat_2));

DA = Aana-Afdj;

DA
printf("merr(DA) %f\n",merr(DA));
tol = 1e-7;
printf("test OK %d\n",merr(DA)<tol);

Aana
Afdj
