## usage: [x3d,ico3d]= extrude (x2d,ico2d,nlay,dx,nodinc)
## extrude a 2D mesh in the z direction. 
## Input args:
##   x2d,ico2d: the 2D mesh
##   nlay number of element layers to extrude 
##   dx increment of z-coordinate between layers (default=1)
##   nodinc = increment in node number between layers (defualt = number
##                          of nodes in the 2D mesh)
## Output args:
##   x3d,ico3d: the resulting 3D mesh

## $Id: extrude.m,v 1.1 2002/10/09 13:38:44 mstorti Exp $
##__INSERT_LICENSE__

function [x3d,ico3d]= extrude (x2d,ico2d,nlay,dx,nodinc)

  if nargin<5
    nodinc=rows(x2d);
  endif
  if nargin<4
    dx=1;
  endif
  if nargin<3
    nlay=1;
  endif

  ## if `dx' is a vector, then they are the new (z) coordinates
  if length(dx) < nlay+1
    dx = dx*(0:nlay)';
  endif

  [nnod,ndim]=size(x2d);
  x3d=zeros((nlay+1)*nnod,ndim+1);
  for k=1:nlay+1
    x3d((k-1)*nnod+(1:nnod),:)=[x2d,dx(k)*ones(nnod,1)];
  endfor

  [nelem,nel]=size(ico2d);
  ico3d=zeros(nelem*nlay,2*nel);
  for k=1:nlay
    ico3d((k-1)*nelem+(1:nelem),:) = [ico2d+(k-1)*nodinc,ico2d+k*nodinc];
  endfor
endfunction
