###key proc.m
### $Id: $

source("mat.output");
Aana = full(Mat_2);
Afdj = full(Mat_3);

all(size(Aana)==[16,16]) || error("bad size");

function ix = indices(nodes,ndof);
  ix = [];
  for k=nodes';
    ix = [ix;(k-1)*ndof+(1:ndof)'];
  endfor
endfunction 

function Am = mirrorv1(A,dim);
  if dim==1;
    Am = A;
    ndof = 4;
    n = rows(A);
    rem(n,ndof)==0 || error("bad n");
    nnod = n/ndof;
    indx = (0:nnod-1)'*ndof+3;
    Am(indx,:) *= -1;
  else
    Am = mirrorv(A',1)';
  endif
endfunction 

function Am = mirrorv(A,varargin);
  Am = A;
  for k=1:length(varargin);
    Am = mirrorv1(Am,varargin{k});
  endfor
endfunction 

function Acond = condens(A);

  if 1
    i1 = indices([1,2]',4);
    i2 = indices([3,4]',4);
    A11 = A(i1,i1);
    A12 = A(i1,i2);
    A21 = A(i2,i1);
    A22 = A(i2,i2);
    A12 = mirrorv(A12,2);
    A21 = mirrorv(A21,1);
    A22 = mirrorv(A22,1,2);
    tol = 1e-6;
    merr(A11-A22)<tol || error("A11!=A22");
    merr(A12-A21)<tol || error("A12!=A21");
    merr(A12-0.5*A11)<tol || error("A12!=0.5*A11");
    Acond = A11+A12+A21+A22;
  else
    Acond = A;
  endif
  
endfunction 

# Aana = condens(Aana);
# Afdj = condens(Afdj);

indx = vec(reshape(1:16,4,4)');
DA = Aana-Afdj;
DA = DA(indx,indx);

DA
printf("merr(DA) %g\n",merr(DA));
Aana
Afdj
