## Copyright (C) 2002 Mario A. Storti
##
## This file is part of Octave.
##__INSERT_LICENSE__
## $Id: proc.m,v 1.1 2002/11/10 03:45:02 mstorti Exp $

##usage: octave> proc

## Author: Mario Storti
## Keywords: DNS, LES, turbulence, periodic flow, plane Poiseuille flow,
## stability

source("data.m.tmp");
load stabi.some_rslt.tmp
U=stabi;
clear stabi

rem(rows(U),Ny+1) == 0 || error("Not correct number of entries");
nstepc = rows(U)/(Ny+1);
printf("%d computed steps\n",nstepc);

u=reshape(U(:,2),Ny+1,nstepc);
v=reshape(U(:,3),Ny+1,nstepc);
w=reshape(U(:,4),Ny+1,nstepc);
p=reshape(U(:,5),Ny+1,nstepc);
