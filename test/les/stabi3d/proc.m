## Copyright (C) 2002 Mario A. Storti
##
## This file is part of Octave.
##__INSERT_LICENSE__
## $Id: proc.m,v 1.7 2002/11/25 20:32:34 mstorti Exp $

##usage: octave> proc

## Author: Mario Storti
## Keywords: DNS, LES, turbulence, periodic flow, plane Poiseuille flow,
## stability

source("data.m.tmp");
if 1
  load stabi.some_rslt.tmp
  U=stabi;
  clear stabi

  load stabi.nod.tmp
  x = stabi;
  clear stabi;

  some = aload("stabi.some.tmp");
  iy = some(1:Ny+1);
  iz=some;
  iz(1:Ny+1)=[];

  y = x(iy,2);
  z = x(iz,3);
  rem(Ny,2)==0 || warning("Sould use Ny even!!");
  iz = [round(Ny/2)+1 Ny+1+(1:Nz)]';
  clear x

endif

if 0
  load stabi.some_rslt4.tmp
  UU = stabi;
  clear stabi;
  U=[UU;U];
endif

## nsome=Ny+1;
nsome=(Ny+Nz+1);
rem(rows(U),nsome) == 0 || error("Not correct number of entries");
nstepc = rows(U)/nsome;

u=reshape(U(:,2),nsome,nstepc);
v=reshape(U(:,3),nsome,nstepc);
w=reshape(U(:,4),nsome,nstepc);
p=reshape(U(:,5),nsome,nstepc);

if 0
  indx=(121:125)';
  u(:,indx)=[];
  v(:,indx)=[];
  w(:,indx)=[];
  p(:,indx)=[];
  nstepc = columns(u);
endif
printf("%d computed steps\n",nstepc);
