%% usage: C= rowprod (A,B,m,n,p)
%%
%% perform matrix products of the form C_k = A_k * B_k where A_k,B_k,C_k
%% are matrices of sizes mxp, pxn, mxn formed with the k-th rows of
%% A,B,C (fortran indexing assumed). If n=1, then m,n and p are
%% optional. In that case you need, off course
%% columns(A)=columns(B)*columns(C), and rows(A)=rows(x)=rows(y)
function C= rowprod (A,B,m,n,p)

%%$Id: rowprod.m,v 1.1.1.1 2000/12/28 12:54:43 mstorti Exp $
  
  [nA,mp]=size(A);
  [nB,pn]=size(B);
  
  if nargin==2
    n=1;
    p=pn;
    m=mp/p;
  end%if
  mn=m*n;

  if mp~=m*p | pn~=p*n
    error('matrix dimensions don''t verify dimension rules');
  end%if
%  n=1;
  
  if nB~=nA
    error('the nomber of rows for A,B,C should be the same')
  end%if
  N=nA;
  C=zeros(N,m*n);

  for k=1:m
    for l=1:n
      sum=zeros(N,1);
      for q=1:p
%	k,l,q
% 	(q-1)*m+k
% 	(l-1)*p+q
% 	(l-1)*m+k
	sum=sum+A(:,(q-1)*m+k).*B(:,(l-1)*p+q);
      end%for
      C(:,(l-1)*m+k)=sum;
    end%for
  end%for
      
