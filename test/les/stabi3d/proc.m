## Copyright (C) 2002 Mario A. Storti
##
## This file is part of Octave.
##__INSERT_LICENSE__
## $Id: proc.m,v 1.2 2002/11/10 14:26:43 mstorti Exp $

##usage: octave> proc

## Author: Mario Storti
## Keywords: DNS, LES, turbulence, periodic flow, plane Poiseuille flow,
## stability

source("data.m.tmp");
load stabi.some_rslt.tmp
U=stabi;
clear stabi

nsome=Ny+1;
nsome=21;
rem(rows(U),nsome) == 0 || error("Not correct number of entries");
nstepc = rows(U)/nsome;
printf("%d computed steps\n",nstepc);

u=reshape(U(:,2),nsome,nstepc);
v=reshape(U(:,3),nsome,nstepc);
w=reshape(U(:,4),nsome,nstepc);
p=reshape(U(:,5),nsome,nstepc);
