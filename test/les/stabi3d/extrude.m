%usage: 
function [ico3d,x3d]= extrude (x2d,ico2d,nlay,dx,nodinc)

  if nargin<5
    nodinc=rows(x2d);
  endif
  if nargin<4
    dx=1;
  endif
  if nargin<3
    nlay=1;
  endif

  [nnod,ndim]=size(x2d);
  x3d=zeros((nlay+1)*nnod,ndim+1);
  for k=1:nlay+1
    x3d((k-1)*nnod+(1:nnod),:)=[x2d,(k-1)*dx*ones(nnod,1)];
  endfor

  [nelem,nel]=size(ico2d);
  ico3d=zeros(nelem*nlay,2*nel);
  for k=1:nlay
    ico3d((k-1)*nelem+(1:nelem),:) = [ico2d+(k-1)*nodinc,ico2d+k*nodinc];
  endfor
endfunction
