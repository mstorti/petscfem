## Copyright (C) 2002 Mario A. Storti
##
## This file is part of Octave.
##__INSERT_LICENSE__
## $Id: tryme4.m,v 1.2 2003/03/04 21:19:57 mstorti Exp $

## Tries to determine if from a set of M points we can find
## the N<M points such that the stencil of these N points have
## the best cloud quality. 

## Author: Mario Storti
## Keywords: cloud, least squares
m=2;
w = zhomo([0 1 0 1],2*m+1,2*m+1);
[X,icone] = pfcm2fem(w);
X = X*2*m;

x=X(:,1); y=X(:,2);
N=rows(X);

A=[ones(N,1) x x.*x  y y.*x y.*x.*x y.*y y.*y.*x y.*y.*x.*x];
A = leftscal(1./l2(A),A);
nA = (1:N)';
indx = 2*(2*m+1)+m+1;
AA=[];
nAA=[];

while 1
  AA=[AA;
      A(indx,:)];
  A(indx,:)=[];
  nAA=[nAA;
       nA(indx)];
  nA(indx)=[];

  if rows(A)==0; break; endif
  H = A*AA';
  [bid,indx] = min(max(H'));
endwhile
